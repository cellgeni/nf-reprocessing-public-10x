#!/usr/bin/env python3
"""
Collect all per-run FASTQ files belonging to ONE sample and rename them into a
single Cell Ranger-style FASTQ set, giving each run its own lane.

Placement: run infer_10x_run.py once per run, then run this on all of that
sample's run FASTQs. The result is SAMPLE_S1_L001_*, SAMPLE_S1_L002_*, ...
ready for STARsolo / Cell Ranger, where every run becomes a distinct lane
(L00n) under one shared sample name and S-index.

Cell Ranger naming: <sample>_S<index>_L<lane:03d>_<readtype>_001.fastq.gz

Read-length policy (compared per read type, across runs):
  * R1 / R2 / R3 (barcode/UMI/biological/genomic reads): lengths MUST agree
    across runs -- a mismatch is fatal.
  * I1 / I2 (sample/cell index reads): a length mismatch is NON-fatal; that
    index read type is dropped from the output with a warning.

How runs are grouped from a flat file list: each input is split into a read-role
token (R1/R2/R3/R4/I1/I2) and the surrounding "stem". A run is identified by
(containing directory, stem), so this works whether the per-run outputs sit in
one flat directory with distinct prefixes, or in per-run subdirectories that
share identical basenames. Use --dry-run first to confirm the detected layout.
"""

from __future__ import annotations

import argparse
import bz2
import gzip
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

ROLE_RE = re.compile(r"(?:^|[_.\-])(R[1-4]|I[12])(?:[_.\-]|$)", re.IGNORECASE)
SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
LANE_WIDTH = 3  # Cell Ranger style: L001


class RenameError(RuntimeError):
    pass


@dataclass
class FileRec:
    path: Path           # original path as given
    real: Path           # resolved path
    basename: str
    role: str            # R1/R2/R3/R4/I1/I2 (uppercased)
    stem: str            # basename with the role token (and trailing chunk) removed
    parent: str          # resolved parent directory, used to separate runs
    length: int          # modal sampled read length


# ---------------------------------------------------------------------------
# Small FASTQ / IO helpers
# ---------------------------------------------------------------------------


def log(level: str, msg: str) -> None:
    print(f"[{level:<5}] {msg}", file=sys.stderr)


def compression_of(path: Path) -> str:
    with open(path, "rb") as fh:
        magic = fh.read(3)
    if magic[:2] == b"\x1f\x8b":
        return "gz"
    if magic[:3] == b"BZh":
        return "bz2"
    return "none"


def open_text_auto(path: Path):
    comp = compression_of(path)
    if comp == "gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    if comp == "bz2":
        return bz2.open(path, "rt", encoding="utf-8", errors="strict")
    return open(path, "rt", encoding="utf-8", errors="strict")


def modal_read_length(path: Path, sample_records: int) -> int:
    """Return the most common sequence length over the first `sample_records`
    reads. Sampling a prefix keeps this cheap regardless of file size."""
    counts: Counter[int] = Counter()
    n = 0
    with open_text_auto(path) as fh:
        while n < sample_records:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()
            if not (s and p and q):
                raise RenameError(f"{path}: truncated FASTQ record at record {n + 1}")
            if n == 0 and not h.startswith("@"):
                raise RenameError(f"{path}: does not look like FASTQ (first line is not a header)")
            counts[len(s.rstrip("\r\n"))] += 1
            n += 1
    if not counts:
        raise RenameError(f"{path}: no FASTQ records found")
    return counts.most_common(1)[0][0]


def split_role(basename: str) -> Tuple[Optional[str], str]:
    """Find the rightmost read-role token; return (role, stem-before-token)."""
    matches = list(ROLE_RE.finditer(basename))
    if not matches:
        return None, ""
    m = matches[-1]
    return m.group(1).upper(), basename[: m.start()]


def natural_key(value: str) -> List[object]:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", value)]


def is_biological(role: str) -> bool:
    # R1/R2/R3/R4 -> fatal on length mismatch; I1/I2 -> skippable.
    return role.startswith("R")


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------


def scan_inputs(paths: Sequence[Path], sample_records: int) -> List[FileRec]:
    recs: List[FileRec] = []
    seen: Dict[Path, Path] = {}
    for path in paths:
        real = path.resolve()
        if real in seen:
            raise RenameError(f"Duplicate input FASTQ after resolving links: {path} == {seen[real]}")
        seen[real] = path
        if compression_of(path) != "gz":
            raise RenameError(
                f"{path}: input must be gzip-compressed (Cell Ranger requires real .fastq.gz). "
                f"gzip it first, or extend this script to recompress on copy."
            )
        role, stem = split_role(path.name)
        if role is None:
            raise RenameError(
                f"{path}: could not find a read-role token (R1/R2/R3/R4/I1/I2) in the filename."
            )
        recs.append(
            FileRec(
                path=path,
                real=real,
                basename=path.name,
                role=role,
                stem=stem,
                parent=str(real.parent),
                length=modal_read_length(path, sample_records),
            )
        )
    return recs


@dataclass
class Assignment:
    rec: FileRec
    lane: int
    run_label: str
    output_name: str


def build_plan(recs: List[FileRec], args: argparse.Namespace) -> Tuple[List[Assignment], List[str], List[str]]:
    warnings: List[str] = []

    # Group files into runs by (directory, stem).
    runs: "defaultdict[Tuple[str, str], List[FileRec]]" = defaultdict(list)
    for r in recs:
        runs[(r.parent, r.stem)].append(r)

    # One file per (run, role).
    for key, members in runs.items():
        by_role: Dict[str, List[FileRec]] = defaultdict(list)
        for r in members:
            by_role[r.role].append(r)
        for role, files in by_role.items():
            if len(files) > 1:
                names = ", ".join(f.basename for f in files)
                raise RenameError(f"Run {key[1] or Path(key[0]).name!r} has multiple {role} files: {names}")

    # Deterministic lane order; build human-friendly, unique run labels.
    ordered_keys = sorted(runs, key=lambda k: (natural_key(k[1]), natural_key(k[0])))
    raw_labels = [k[1] if k[1] else Path(k[0]).name for k in ordered_keys]
    label_counts = Counter(raw_labels)
    labels: Dict[Tuple[str, str], str] = {}
    for key, raw in zip(ordered_keys, raw_labels):
        labels[key] = raw if label_counts[raw] == 1 else f"{raw} ({Path(key[0]).name})"

    lanes: Dict[Tuple[str, str], int] = {k: args.start_lane + i for i, k in enumerate(ordered_keys)}

    # Sanity warning: biological read sets that differ between runs.
    bio_sets = {labels[k]: tuple(sorted({r.role for r in runs[k] if is_biological(r.role)})) for k in ordered_keys}
    if len(set(bio_sets.values())) > 1:
        detail = "; ".join(f"{lab}:{'+'.join(s) or 'none'}" for lab, s in bio_sets.items())
        warnings.append(f"Runs do not share the same set of biological reads: {detail}")

    # Per-role length comparison across runs.
    by_role_lengths: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    for key in ordered_keys:
        for r in runs[key]:
            by_role_lengths[r.role].append((labels[key], r.length))

    skipped_roles: List[str] = []
    for role in sorted(by_role_lengths):
        items = by_role_lengths[role]
        if len({length for _, length in items}) <= 1:
            continue
        detail = ", ".join(f"{lab}={length}" for lab, length in items)
        if is_biological(role):
            raise RenameError(f"{role} read length differs across runs: {detail}")
        skipped_roles.append(role)
        warnings.append(f"{role} read length differs across runs ({detail}); skipping all {role} files")

    # Build assignments for the kept roles.
    assignments: List[Assignment] = []
    for key in ordered_keys:
        lane = lanes[key]
        for r in runs[key]:
            if r.role in skipped_roles:
                continue
            out = f"{args.sample_id}_S{args.sample_index}_L{lane:0{LANE_WIDTH}d}_{r.role}_001.fastq.gz"
            assignments.append(Assignment(rec=r, lane=lane, run_label=labels[key], output_name=out))

    out_names = [a.output_name for a in assignments]
    dups = [n for n, c in Counter(out_names).items() if c > 1]
    if dups:
        raise RenameError(f"Duplicate output names in plan (input layout is ambiguous): {dups}")

    return assignments, skipped_roles, warnings


def place(src: Path, dest: Path, action: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if action == "symlink":
        os.symlink(os.path.abspath(src), dest)
    elif action == "hardlink":
        try:
            os.link(os.path.abspath(src), dest)
        except OSError:
            shutil.copy2(src, dest)
    else:  # copy
        shutil.copy2(src, dest)


def write_manifest(assignments: List[Assignment], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        out.write("run_label\tlane\tread_type\tread_length\tsource\toutput_name\n")
        for a in assignments:
            out.write(
                f"{a.run_label}\t{a.lane}\t{a.rec.role}\t{a.rec.length}\t{a.rec.path}\t{a.output_name}\n"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Aggregate one sample's per-run FASTQs and rename to Cell Ranger convention (one lane per run)."
    )
    p.add_argument("--fastqs", nargs="+", type=Path, required=True, help="All FASTQ.GZ files for one sample, across runs")
    p.add_argument("--sample-id", required=True, help="Sample name used in the Cell Ranger prefix")
    p.add_argument("--sample-index", type=int, default=1, help="S-index (default 1)")
    p.add_argument("--outdir", type=Path, required=True, help="Where renamed FASTQs are written")
    p.add_argument("--action", choices=("symlink", "hardlink", "copy"), default="hardlink",
                   help="How to materialise outputs (default hardlink)")
    p.add_argument("--start-lane", type=int, default=1, help="Lane number for the first run (default 1)")
    p.add_argument("--sample-records", type=int, default=1000, help="Reads sampled per file to get its length")
    p.add_argument("--manifest", type=Path, default=None, help="Optional TSV mapping of source -> output")
    p.add_argument("--allow-unsafe-name", action="store_true", help="Permit sample names outside [A-Za-z0-9_.-]")
    p.add_argument("--dry-run", action="store_true", help="Plan and report, but do not write any files")
    args = p.parse_args(argv)

    for f in args.fastqs:
        if not f.exists():
            p.error(f"FASTQ not found: {f}")
    if not args.allow_unsafe_name and not SAFE_NAME_RE.match(args.sample_id):
        p.error(f"--sample-id {args.sample_id!r} has characters unsafe for Cell Ranger names; "
                f"use --allow-unsafe-name to override")
    if args.start_lane < 1:
        p.error("--start-lane must be >= 1")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        recs = scan_inputs(args.fastqs, args.sample_records)
        assignments, skipped, warnings = build_plan(recs, args)

        n_runs = len({(a.run_label, a.lane) for a in assignments})
        log("INFO", f"Sample {args.sample_id}: {len(recs)} files across {n_runs} run(s)")
        for lane, label in sorted({(a.lane, a.run_label) for a in assignments}):
            roles = sorted(a.rec.role for a in assignments if a.lane == lane)
            log("INFO", f"  L{lane:0{LANE_WIDTH}d}  <- {label}  [{', '.join(roles)}]")
        for w in warnings:
            log("WARN", w)
        if skipped:
            log("WARN", f"Skipped read type(s) due to length mismatch: {', '.join(skipped)}")

        for a in assignments:
            log("INFO", f"  {a.rec.basename}  ->  {a.output_name}  (len={a.rec.length})")
            if not args.dry_run:
                place(a.rec.path, args.outdir / a.output_name, args.action)

        if args.manifest and not args.dry_run:
            write_manifest(assignments, args.manifest)

        log("INFO", f"{'Planned' if args.dry_run else 'Wrote'} {len(assignments)} output FASTQ(s)"
                    + (" (dry-run)" if args.dry_run else f" to {args.outdir}"))
        return 0
    except RenameError as exc:
        log("ERROR", str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())