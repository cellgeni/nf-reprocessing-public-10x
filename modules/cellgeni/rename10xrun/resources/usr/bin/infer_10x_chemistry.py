#!/usr/bin/env python3
"""
Infer 10x FASTQ read roles and normalize filenames for Cell Ranger/STARsolo.

This tool is intentionally conservative. It emits Cell Ranger-style FASTQs only
when barcode whitelist evidence and read-length geometry identify a unique 10x
layout. Ambiguous or mixed-modality inputs fail closed unless --prefer-gex is
used to select the GEX library from a mixed Multiome sample. ATAC-only inputs
are identified but are not emitted unless --allow-atac-output is explicit.
"""

from __future__ import annotations

import argparse
import bz2
import collections
import gzip
import json
import os
import re
import shutil
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

DNA_RE = re.compile(r"^[ACGTNacgtn]+$")
LANE_RE = re.compile(r"(?:^|[_\-.])L00([1-8])(?:[_\-.]|$)")
READ_RE = re.compile(r"(?:^|[_\-.])(R1|R2|R3|I1|I2|_1|_2|_3|1|2|3)(?:[_\-.]|$)", re.IGNORECASE)
SAMPLE_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


@dataclass(frozen=True)
class WhitelistDef:
    id: str
    label: str
    cb_len: int
    filenames: Tuple[str, ...]
    modality: str
    umi_len: Optional[int] = None
    usable_for_starsolo: bool = True


WHITELISTS: Tuple[WhitelistDef, ...] = (
    WhitelistDef(
        id="10x_3p_v1",
        label="Chromium Single Cell 3prime v1",
        cb_len=14,
        filenames=("737K-april-2014_rc.txt", "737K-april-2014_rc.txt.gz"),
        modality="gex",
        umi_len=None,
    ),
    WhitelistDef(
        id="10x_3p_v2_or_5p_v1_v2",
        label="Chromium 3prime v2 or 5prime v1/v2",
        cb_len=16,
        filenames=("737K-august-2016.txt", "737K-august-2016.txt.gz"),
        modality="gex",
        umi_len=10,
    ),
    WhitelistDef(
        id="10x_3p_v3_family",
        label="Chromium 3prime v3/v3.1/LT/HT",
        cb_len=16,
        filenames=("3M-february-2018.txt", "3M-february-2018.txt.gz"),
        modality="gex",
        umi_len=12,
    ),
    WhitelistDef(
        id="10x_3p_v4_gemx",
        label="GEM-X Universal 3prime v4",
        cb_len=16,
        filenames=("3M-3pgex-may-2023.txt", "3M-3pgex-may-2023.txt.gz"),
        modality="gex",
        umi_len=12,
    ),
    WhitelistDef(
        id="10x_5p_v3_gemx",
        label="GEM-X Universal 5prime v3",
        cb_len=16,
        filenames=("3M-5pgex-jan-2023.txt", "3M-5pgex-jan-2023.txt.gz"),
        modality="gex",
        umi_len=12,
    ),
    WhitelistDef(
        id="10x_multiome_arc_gex",
        label="Chromium Single Cell Multiome Gene Expression",
        cb_len=16,
        filenames=(
            "gex_737K-arc-v1.txt",
            "gex_737K-arc-v1.txt.gz",
            "737K-arc-v1.txt",
            "737K-arc-v1.txt.gz",
        ),
        modality="gex",
        umi_len=12,
    ),
    WhitelistDef(
        id="10x_multiome_arc_atac",
        label="Chromium Single Cell Multiome ATAC",
        cb_len=16,
        filenames=(
            "atac_737K-arc-v1.txt",
            "atac_737K-arc-v1.txt.gz",
            "737K-arc-v1.txt",
            "737K-arc-v1.txt.gz",
        ),
        modality="atac",
        usable_for_starsolo=False,
    ),
    WhitelistDef(
        id="10x_atac_v1_v1.1_v2",
        label="Chromium Single Cell ATAC v1/v1.1/v2",
        cb_len=16,
        filenames=("737K-cratac-v1.txt", "737K-cratac-v1.txt.gz"),
        modality="atac",
        usable_for_starsolo=False,
    ),
)


@dataclass
class MatchResult:
    whitelist_id: str
    label: str
    modality: str
    cb_len: int
    umi_len: Optional[int]
    offset: int
    hits: int
    sampled: int
    fraction: float


@dataclass
class FastqStats:
    path: str
    basename: str
    sampled_records: int
    length_counts: Dict[int, int]
    common_length: int
    common_fraction: float
    min_length: int
    max_length: int
    existing_read_type: Optional[str] = None
    existing_lane: Optional[str] = None
    matches: List[MatchResult] = field(default_factory=list)

    def to_json(self) -> Dict[str, object]:
        data = asdict(self)
        data["length_counts"] = {str(k): v for k, v in sorted(self.length_counts.items())}
        return data


@dataclass
class OutputEntry:
    read_type: str
    lane: str
    output_name: str
    sources: List[str]
    synthetic: bool = False
    umi_len: Optional[int] = None


@dataclass
class InferencePlan:
    sample: str
    chemistry_id: str
    chemistry_label: str
    modality: str
    usable_for_starsolo: bool
    confidence: str
    warnings: List[str]
    entries: List[OutputEntry]
    excluded_files: List[str]
    stats: List[FastqStats]
    checks: Dict[str, object] = field(default_factory=dict)


class TenxInferenceError(RuntimeError):
    pass


class WhitelistStore:
    def __init__(self, whitelist_dir: Path):
        self.whitelist_dir = whitelist_dir
        self._sets: Dict[str, Optional[set[str]]] = {}
        self._paths: Dict[str, Optional[Path]] = {}
        for spec in WHITELISTS:
            self._paths[spec.id] = self._find_file(spec)

    def _find_file(self, spec: WhitelistDef) -> Optional[Path]:
        for name in spec.filenames:
            p = self.whitelist_dir / name
            if p.exists():
                return p
        return None

    def available_specs(self) -> List[WhitelistDef]:
        return [spec for spec in WHITELISTS if self._paths.get(spec.id) is not None]

    def get(self, spec: WhitelistDef) -> Optional[set[str]]:
        if spec.id in self._sets:
            return self._sets[spec.id]
        path = self._paths.get(spec.id)
        if path is None:
            self._sets[spec.id] = None
            return None
        values: set[str] = set()
        with open_text_auto(path) as handle:
            for line in handle:
                bc = line.strip().split()[0] if line.strip() else ""
                if len(bc) >= spec.cb_len:
                    values.add(bc[: spec.cb_len].upper())
        self._sets[spec.id] = values
        return values


def compression_of(path: Path) -> str:
    with open(path, "rb") as handle:
        magic = handle.read(3)
    if magic[:2] == b"\x1f\x8b":
        return "gz"
    if magic[:3] == b"BZh":
        return "bz2"
    return "none"


def open_text_auto(path: Path):
    comp = compression_of(path)
    if comp == "gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if comp == "bz2":
        return bz2.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def open_bytes_auto(path: Path):
    comp = compression_of(path)
    if comp == "gz":
        return gzip.open(path, "rb")
    if comp == "bz2":
        return bz2.open(path, "rb")
    return open(path, "rb")


def iter_fastq_records(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with open_text_auto(path) as handle:
        recno = 0
        while True:
            h = handle.readline()
            if not h:
                break
            s = handle.readline()
            p = handle.readline()
            q = handle.readline()
            recno += 1
            if not (s and p and q):
                raise TenxInferenceError(f"{path}: truncated FASTQ record at record {recno}")
            if not h.startswith("@"):
                raise TenxInferenceError(f"{path}: FASTQ header at record {recno} does not start with @")
            if not p.startswith("+"):
                raise TenxInferenceError(f"{path}: FASTQ plus line at record {recno} does not start with +")
            seq = s.rstrip("\r\n")
            qual = q.rstrip("\r\n")
            if len(seq) != len(qual):
                raise TenxInferenceError(f"{path}: sequence/quality length mismatch at record {recno}")
            if not DNA_RE.match(seq):
                raise TenxInferenceError(f"{path}: non-DNA sequence characters at record {recno}")
            yield h.rstrip("\r\n"), seq, p.rstrip("\r\n"), qual


def normalized_read_id(header: str) -> str:
    x = header[1:] if header.startswith("@") else header
    x = x.split()[0]
    x = re.sub(r"([/.][123])$", "", x)
    return x


def natural_key(path: str) -> Tuple[object, ...]:
    parts = re.split(r"(\d+)", Path(path).name)
    return tuple(int(p) if p.isdigit() else p.lower() for p in parts)


def parse_existing_read_type(name: str) -> Optional[str]:
    # Prefer explicit Cell Ranger-style read tokens. R3 is common in public
    # Multiome/ATAC exports where source R2 is the barcode index and source R3
    # is the second genomic read that must become Cell Ranger R2.
    m = re.search(r"(?:^|[_\-.])(R1|R2|R3|I1|I2)(?:[_\-.]|$)", name, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    # Common SRA dump suffixes: _1/_2/_3. Treat as hints only.
    m = re.search(r"(?:^|[_\-.])([123])(?:\.f(?:ast)?q(?:\.gz|\.bz2)?|\.fq(?:\.gz|\.bz2)?|[_\-.]|$)", name, re.IGNORECASE)
    if m:
        return "R" + m.group(1)
    return None


def parse_lane(name: str) -> Optional[str]:
    m = LANE_RE.search(name)
    if not m:
        return None
    return "L00" + m.group(1)


def collect_stats(paths: Sequence[Path], sample_records: int) -> List[FastqStats]:
    stats: List[FastqStats] = []
    for path in paths:
        counts: collections.Counter[int] = collections.Counter()
        n = 0
        for _h, seq, _p, _q in iter_fastq_records(path):
            counts[len(seq)] += 1
            n += 1
            if n >= sample_records:
                break
        if n == 0:
            raise TenxInferenceError(f"{path}: no FASTQ records sampled")
        common_length, common_n = counts.most_common(1)[0]
        stats.append(
            FastqStats(
                path=str(path),
                basename=path.name,
                sampled_records=n,
                length_counts=dict(counts),
                common_length=common_length,
                common_fraction=common_n / n,
                min_length=min(counts),
                max_length=max(counts),
                existing_read_type=parse_existing_read_type(path.name),
                existing_lane=parse_lane(path.name),
            )
        )
    return stats


def sample_sequences(path: Path, limit: int) -> List[str]:
    seqs: List[str] = []
    for _h, seq, _p, _q in iter_fastq_records(path):
        seqs.append(seq.upper())
        if len(seqs) >= limit:
            break
    return seqs


def match_whitelists(stats: List[FastqStats], wl_store: WhitelistStore, args: argparse.Namespace) -> None:
    available = wl_store.available_specs()
    missing = [spec for spec in WHITELISTS if wl_store._paths.get(spec.id) is None]
    if missing and not args.allow_missing_whitelists:
        # Missing RNA/GEX whitelists make exact chemistry inference incomplete and
        # are fatal. Missing ATAC-only whitelists are reported but not fatal because
        # Multiome ATAC is also detected using the shared ARC whitelist and ATAC-like
        # geometry is rejected rather than processed as GEX.
        missing_required = [spec for spec in missing if spec.modality == "gex"]
        detail = "; ".join(f"{spec.id}: one of {','.join(spec.filenames)}" for spec in missing)
        if missing_required:
            raise TenxInferenceError(
                f"Missing required whitelist(s) under {wl_store.whitelist_dir}: {detail}. "
                "Use --allow-missing-whitelists only for tests or explicitly limited deployments."
            )
        print(f"WARNING: Missing optional ATAC whitelist(s) under {wl_store.whitelist_dir}: {detail}", file=sys.stderr)
    if not available:
        raise TenxInferenceError(
            f"No known whitelist files found under {wl_store.whitelist_dir}. Expected one of: "
            + ", ".join(sorted({x for spec in WHITELISTS for x in spec.filenames}))
        )
    for st in stats:
        seqs = sample_sequences(Path(st.path), args.sample_records)
        for spec in available:
            wl = wl_store.get(spec)
            if not wl:
                continue
            offsets = [0]
            # Multiome ATAC i5 can be 24 nt with an 8 nt dark/spacer segment. Try
            # both offset 0 and offset 8 so the report is explicit about evidence.
            if spec.modality == "atac" and st.common_length >= 24:
                offsets = [0, 8]
            for offset in offsets:
                if st.common_length < offset + spec.cb_len:
                    continue
                hits = 0
                for seq in seqs:
                    if len(seq) >= offset + spec.cb_len and seq[offset : offset + spec.cb_len] in wl:
                        hits += 1
                frac = hits / max(1, len(seqs))
                if hits >= args.min_whitelist_hits and frac >= args.min_whitelist_fraction:
                    st.matches.append(
                        MatchResult(
                            whitelist_id=spec.id,
                            label=spec.label,
                            modality=spec.modality,
                            cb_len=spec.cb_len,
                            umi_len=spec.umi_len,
                            offset=offset,
                            hits=hits,
                            sampled=len(seqs),
                            fraction=frac,
                        )
                    )
        st.matches.sort(key=lambda m: (m.fraction, m.hits), reverse=True)


def best_match(st: FastqStats, modality: Optional[str] = None) -> Optional[MatchResult]:
    matches = [m for m in st.matches if modality is None or m.modality == modality]
    return matches[0] if matches else None


def require_constant_length(st: FastqStats, role: str, min_fraction: float = 0.98) -> None:
    if st.common_fraction < min_fraction:
        raise TenxInferenceError(
            f"{st.basename}: {role} has variable read lengths "
            f"{st.length_counts}; this looks trimmed or mixed and is unsafe to rename"
        )


def is_index_like(st: FastqStats) -> bool:
    return st.common_length in {6, 7, 8, 9, 10, 14, 16, 24} and st.common_fraction >= 0.95


def is_short_rna_barcode(st: FastqStats, m: MatchResult) -> bool:
    if m.whitelist_id == "10x_3p_v1":
        # Either separate v1 CB read (14 nt) or already merged CB+UMI read.
        return st.common_length in {14, 19, 24} or 18 <= st.common_length <= 26
    expected = m.cb_len + (m.umi_len or 0)
    return expected <= st.common_length <= max(expected + 4, 32)


def classify_gex(stats: List[FastqStats], args: argparse.Namespace, allow_excluded: bool) -> InferencePlan:
    warnings: List[str] = []
    gex_candidates = [st for st in stats if has_gex_layout(st)]
    if not gex_candidates:
        raise TenxInferenceError("No GEX-compatible barcode read matched any 10x RNA whitelist")

    # Pick a chemistry by the best matching whitelist. All barcode reads must be consistent.
    candidates_by_id: Dict[str, List[FastqStats]] = collections.defaultdict(list)
    for st in gex_candidates:
        m = best_match(st, "gex")
        if m is None:
            continue
        if has_gex_layout(st):
            candidates_by_id[m.whitelist_id].append(st)
    if not candidates_by_id:
        raise TenxInferenceError("Whitelist matches were present but not on reads with RNA barcode geometry")

    best_id = max(candidates_by_id, key=lambda k: sum(best_match(s, "gex").hits for s in candidates_by_id[k]))
    r1_files = sorted(candidates_by_id[best_id], key=lambda st: natural_key(st.path))
    first_match = best_match(r1_files[0], "gex")
    assert first_match is not None

    for st in r1_files:
        m = best_match(st, "gex")
        if m is None or m.whitelist_id != best_id:
            raise TenxInferenceError("Inconsistent RNA barcode whitelist matches across candidate R1 files")
        require_constant_length(st, "RNA barcode read")

    # Known ambiguity that does not affect naming or STARsolo CB/UMI positions.
    confidence = "unique"
    if best_id == "10x_3p_v2_or_5p_v1_v2":
        confidence = "layout_only"
        warnings.append("737K-august-2016 cannot distinguish 3prime v2 from 5prime v1/v2 by FASTQ geometry alone")
    if best_id == "10x_3p_v3_family":
        confidence = "layout_only"
        warnings.append("3prime v3, v3.1, LT and HT share the same whitelist/layout family")

    remaining = [st for st in stats if st not in r1_files]
    excluded: List[str] = []

    # v1 can be delivered as separate CB, UMI and biological FASTQs.
    if best_id == "10x_3p_v1" and all(st.common_length == 14 for st in r1_files):
        return classify_v1_separate(stats, r1_files, first_match, args)

    expected_len = first_match.cb_len + (first_match.umi_len or 0)
    umi_len = first_match.umi_len
    if best_id == "10x_3p_v1" and first_match.umi_len is None:
        umi_len = r1_files[0].common_length - first_match.cb_len
        if umi_len not in args.allow_v1_umi_lengths:
            raise TenxInferenceError(
                f"Inferred 3prime v1 merged CB+UMI length {r1_files[0].common_length} gives UMI length {umi_len}, "
                f"not in allowed values {sorted(args.allow_v1_umi_lengths)}"
            )
    else:
        too_short = [st.basename for st in r1_files if st.common_length < expected_len]
        if too_short:
            raise TenxInferenceError(f"RNA barcode reads shorter than CB+UMI ({expected_len}): {too_short}")
        too_long = [st.basename for st in r1_files if expected_len < st.common_length < 50]
        if too_long:
            warnings.append(f"RNA barcode reads include bases beyond CB+UMI ({expected_len}): {too_long}")

    # Biological reads. For standard single-end RNA layout, R2 is the only long read.
    # For mixed Multiome, --prefer-gex selects the longer GEX R2 and excludes ATAC reads.
    long_reads = [st for st in remaining if st.common_length >= args.min_bio_read_length]
    if args.prefer_gex:
        gex_long = [st for st in long_reads if st.common_length >= args.min_gex_r2_length]
        if len(gex_long) >= len(r1_files):
            excluded.extend([st.path for st in long_reads if st not in gex_long])
            long_reads = gex_long
    if len(long_reads) != len(r1_files):
        raise TenxInferenceError(
            f"Expected {len(r1_files)} GEX biological R2 FASTQ(s) for RNA layout, found {len(long_reads)}. "
            "This is ambiguous or mixed modality. Use --prefer-gex only when a Multiome GEX library should be selected."
        )
    r2_files = sorted(long_reads, key=lambda st: natural_key(st.path))
    for st in r2_files:
        require_constant_length(st, "biological read", min_fraction=0.90)

    used_paths = {st.path for st in (r1_files + r2_files)}
    index_files = [st for st in stats if st.path not in used_paths and st.path not in excluded]
    if args.prefer_gex and any(best_match(st, "atac") is not None for st in stats):
        # Mixed Multiome input: keep the 10-cycle GEX indexes and exclude ATAC
        # 8/16/24-cycle index reads so they cannot be passed to STARsolo with GEX.
        gex_index_files = [st for st in index_files if st.common_length == 10]
        excluded.extend([st.path for st in index_files if st not in gex_index_files])
        index_files = gex_index_files
    entries: List[OutputEntry] = []
    lanes = assign_lanes({"R1": r1_files, "R2": r2_files})
    for lane, r1, r2 in zip(lanes, r1_files, r2_files):
        entries.append(OutputEntry("R1", lane, cellranger_name(args.sample, lane, "R1"), [r1.path], umi_len=umi_len))
        entries.append(OutputEntry("R2", lane, cellranger_name(args.sample, lane, "R2"), [r2.path]))

    add_index_entries(entries, index_files, args.sample, lanes, warnings)

    if index_files and len(index_files) not in {0, len(r1_files), 2 * len(r1_files)}:
        leftovers = [st.path for st in index_files]
        if allow_excluded:
            excluded.extend(leftovers)
        else:
            raise TenxInferenceError(f"Unpaired index-like files remain after RNA assignment: {leftovers}")

    fail_if_unassigned(stats, entries, excluded, "RNA/GEX role assignment")

    return InferencePlan(
        sample=args.sample,
        chemistry_id=best_id,
        chemistry_label=first_match.label,
        modality="gex",
        usable_for_starsolo=True,
        confidence=confidence,
        warnings=warnings,
        entries=entries,
        excluded_files=excluded,
        stats=stats,
    )


def classify_v1_separate(
    stats: List[FastqStats], cb_files: List[FastqStats], match: MatchResult, args: argparse.Namespace
) -> InferencePlan:
    warnings: List[str] = []
    n = len(cb_files)
    remaining = [st for st in stats if st not in cb_files]
    umi_files = [st for st in remaining if st.common_length in args.allow_v1_umi_lengths and st.common_fraction >= 0.98]
    bio_files = [st for st in remaining if st.common_length >= args.min_bio_read_length]
    if len(umi_files) != n:
        raise TenxInferenceError(
            f"Detected 3prime v1 14 nt CB read(s), but found {len(umi_files)} UMI FASTQ(s); "
            f"allowed UMI lengths are {sorted(args.allow_v1_umi_lengths)}"
        )
    if len(bio_files) != n:
        raise TenxInferenceError(f"Detected 3prime v1 CB/UMI reads, but found {len(bio_files)} biological FASTQ(s)")
    cb_files = sorted(cb_files, key=lambda st: natural_key(st.path))
    umi_files = sorted(umi_files, key=lambda st: natural_key(st.path))
    bio_files = sorted(bio_files, key=lambda st: natural_key(st.path))
    umi_lengths = {st.common_length for st in umi_files}
    if len(umi_lengths) != 1:
        raise TenxInferenceError(f"3prime v1 UMI reads have inconsistent lengths: {sorted(umi_lengths)}")
    umi_len = next(iter(umi_lengths))
    if umi_len == 5:
        warnings.append("3prime v1 UMI length is 5 nt; this is rare but explicitly allowed")
    entries: List[OutputEntry] = []
    lanes = assign_lanes({"CB": cb_files, "UMI": umi_files, "R2": bio_files})
    for lane, cb, umi, bio in zip(lanes, cb_files, umi_files, bio_files):
        entries.append(
            OutputEntry(
                "R1",
                lane,
                cellranger_name(args.sample, lane, "R1"),
                [cb.path, umi.path],
                synthetic=True,
                umi_len=umi_len,
            )
        )
        entries.append(OutputEntry("R2", lane, cellranger_name(args.sample, lane, "R2"), [bio.path]))
    used_paths = {st.path for st in (cb_files + umi_files + bio_files)}
    index_files = [st for st in stats if st.path not in used_paths]
    add_index_entries(entries, index_files, args.sample, lanes, warnings)
    fail_if_unassigned(stats, entries, [], "3prime v1 role assignment")
    return InferencePlan(
        sample=args.sample,
        chemistry_id="10x_3p_v1",
        chemistry_label="Chromium Single Cell 3prime v1",
        modality="gex",
        usable_for_starsolo=True,
        confidence="unique",
        warnings=warnings,
        entries=entries,
        excluded_files=[],
        stats=stats,
    )


def classify_atac(stats: List[FastqStats], args: argparse.Namespace) -> InferencePlan:
    warnings: List[str] = []
    i2_candidates = [st for st in stats if best_match(st, "atac") is not None and st.common_length in {16, 24}]
    if not i2_candidates:
        raise TenxInferenceError("No ATAC I2 barcode read matched a 10x ATAC/Multiome ATAC whitelist")
    # Pick all I2 candidates with the same best whitelist as the strongest one.
    strongest = max(i2_candidates, key=lambda st: best_match(st, "atac").hits)
    strongest_match = best_match(strongest, "atac")
    assert strongest_match is not None
    i2_files = sorted(
        [st for st in i2_candidates if best_match(st, "atac").whitelist_id == strongest_match.whitelist_id],
        key=lambda st: natural_key(st.path),
    )
    for st in i2_files:
        require_constant_length(st, "ATAC barcode index")
    n = len(i2_files)
    remaining = [st for st in stats if st not in i2_files]
    long_reads = sorted([st for st in remaining if st.common_length >= args.min_atac_genomic_read_length], key=lambda st: natural_key(st.path))
    if len(long_reads) != 2 * n:
        raise TenxInferenceError(f"Expected {2*n} ATAC genomic read FASTQ(s), found {len(long_reads)}")
    # Assign R1/R2. 50/49 is unique for Multiome ATAC. Equal-length ATAC needs existing read labels.
    r1_files: List[FastqStats] = []
    r2_files: List[FastqStats] = []
    per_lane = [long_reads[i : i + 2] for i in range(0, len(long_reads), 2)]
    for pair in per_lane:
        explicit_r1 = [st for st in pair if st.existing_read_type == "R1"]
        explicit_r2 = [st for st in pair if st.existing_read_type == "R2"]
        explicit_r3 = [st for st in pair if st.existing_read_type == "R3"]
        if len(explicit_r1) == 1 and len(explicit_r3) == 1:
            # Public 10x ATAC/Multiome often appears as R1=genomic, R2=barcode, R3=genomic.
            # The genomic R3 must be renamed to Cell Ranger R2.
            r1_files.append(explicit_r1[0])
            r2_files.append(explicit_r3[0])
        elif len(explicit_r1) == 1 and len(explicit_r2) == 1:
            r1_files.append(explicit_r1[0])
            r2_files.append(explicit_r2[0])
        else:
            by_len = sorted(pair, key=lambda st: st.common_length, reverse=True)
            if by_len[0].common_length == by_len[1].common_length:
                raise TenxInferenceError(
                    "ATAC R1/R2 have equal length and no reliable R1/R2 labels; refusing to guess orientation"
                )
            r1_files.append(by_len[0])
            r2_files.append(by_len[1])
    used_paths = {st.path for st in (i2_files + r1_files + r2_files)}
    index_files = [st for st in stats if st.path not in used_paths]
    lanes = assign_lanes({"R1": r1_files, "R2": r2_files, "I2": i2_files})
    entries: List[OutputEntry] = []
    for lane, r1, r2, i2 in zip(lanes, r1_files, r2_files, i2_files):
        entries.append(OutputEntry("R1", lane, cellranger_name(args.sample, lane, "R1"), [r1.path]))
        entries.append(OutputEntry("R2", lane, cellranger_name(args.sample, lane, "R2"), [r2.path]))
        entries.append(OutputEntry("I2", lane, cellranger_name(args.sample, lane, "I2"), [i2.path]))
    add_index_entries(entries, index_files, args.sample, lanes, warnings, force_i1=True)
    fail_if_unassigned(stats, entries, [], "ATAC role assignment")
    return InferencePlan(
        sample=args.sample,
        chemistry_id=strongest_match.whitelist_id,
        chemistry_label=strongest_match.label,
        modality="atac",
        usable_for_starsolo=False,
        confidence="unique",
        warnings=warnings,
        entries=entries,
        excluded_files=[],
        stats=stats,
    )


def assign_lanes(role_files: Dict[str, List[FastqStats]]) -> List[str]:
    counts = {role: len(files) for role, files in role_files.items()}
    non_zero = {n for n in counts.values() if n != 0}
    if len(non_zero) != 1:
        raise TenxInferenceError(f"Read roles do not have matching lane counts: {counts}")
    n = non_zero.pop()
    lane_sets = []
    for files in role_files.values():
        if all(st.existing_lane for st in files):
            lane_sets.append([st.existing_lane for st in files])
    if lane_sets and all(set(x) == set(lane_sets[0]) for x in lane_sets):
        return sorted(lane_sets[0])
    return [f"L{i:03d}" for i in range(1, n + 1)]


def cellranger_name(sample: str, lane: str, read_type: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample).strip("_") or "sample"
    return f"{safe}_S1_{lane}_{read_type}_001.fastq.gz"


def add_index_entries(
    entries: List[OutputEntry],
    index_files: List[FastqStats],
    sample: str,
    lanes: List[str],
    warnings: List[str],
    force_i1: bool = False,
) -> None:
    if not index_files:
        return
    index_files = sorted(index_files, key=lambda st: natural_key(st.path))
    n = len(lanes)
    explicit_i1 = [st for st in index_files if st.existing_read_type == "I1"]
    explicit_i2 = [st for st in index_files if st.existing_read_type == "I2"]
    if force_i1:
        candidates = explicit_i1 or index_files
        if len(candidates) == n:
            for lane, st in zip(lanes, candidates):
                entries.append(OutputEntry("I1", lane, cellranger_name(sample, lane, "I1"), [st.path]))
        elif candidates:
            warnings.append(f"Could not assign ATAC I1 files unambiguously: {[st.basename for st in candidates]}")
        return
    if explicit_i1 and len(explicit_i1) == n:
        for lane, st in zip(lanes, sorted(explicit_i1, key=lambda st: natural_key(st.path))):
            entries.append(OutputEntry("I1", lane, cellranger_name(sample, lane, "I1"), [st.path]))
    if explicit_i2 and len(explicit_i2) == n:
        for lane, st in zip(lanes, sorted(explicit_i2, key=lambda st: natural_key(st.path))):
            entries.append(OutputEntry("I2", lane, cellranger_name(sample, lane, "I2"), [st.path]))
    assigned_paths = {st.path for st in (explicit_i1 + explicit_i2)}
    leftover = [st for st in index_files if st.path not in assigned_paths and is_index_like(st)]
    if not leftover:
        return
    if len(leftover) == n:
        for lane, st in zip(lanes, leftover):
            entries.append(OutputEntry("I1", lane, cellranger_name(sample, lane, "I1"), [st.path]))
        warnings.append("Index read assigned as I1 from length/name order; STARsolo does not use index reads")
    elif len(leftover) == 2 * n:
        warnings.append("Dual index reads assigned as I1/I2 from sorted order because names lacked I1/I2 labels")
        for i, lane in enumerate(lanes):
            st1 = leftover[2 * i]
            st2 = leftover[2 * i + 1]
            entries.append(OutputEntry("I1", lane, cellranger_name(sample, lane, "I1"), [st1.path]))
            entries.append(OutputEntry("I2", lane, cellranger_name(sample, lane, "I2"), [st2.path]))
    else:
        warnings.append(f"Unassigned index-like FASTQs: {[st.basename for st in leftover]}")


def has_gex_layout(st: FastqStats) -> bool:
    m = best_match(st, "gex")
    if m is None:
        return False
    return is_short_rna_barcode(st, m) or st.common_length >= 50


def has_atac_layout(st: FastqStats) -> bool:
    return best_match(st, "atac") is not None and st.common_length in {16, 24}


def entry_source_paths(entries: Sequence[OutputEntry]) -> set[str]:
    return {src for entry in entries for src in entry.sources}


def fail_if_unassigned(stats: List[FastqStats], entries: Sequence[OutputEntry], excluded: Sequence[str], context: str) -> None:
    assigned = entry_source_paths(entries)
    excluded_set = set(excluded)
    unassigned = [st for st in stats if st.path not in assigned and st.path not in excluded_set]
    if unassigned:
        detail = ", ".join(f"{st.basename}({st.common_length}bp,{st.existing_read_type or '?'})" for st in unassigned)
        raise TenxInferenceError(f"Unassigned FASTQ(s) remain after {context}: {detail}. Refusing to guess or silently drop reads.")


def apply_forced_lane(plan: InferencePlan, lane_number: Optional[int]) -> None:
    if lane_number is None:
        return
    existing = {entry.lane for entry in plan.entries}
    if len(existing) > 1:
        raise TenxInferenceError("--lane can only be used when the input resolves to exactly one lane")
    lane = f"L{lane_number:03d}"
    for entry in plan.entries:
        entry.lane = lane
        entry.output_name = cellranger_name(plan.sample, lane, entry.read_type)


def infer_plan(stats: List[FastqStats], args: argparse.Namespace) -> InferencePlan:
    has_gex = any(has_gex_layout(st) for st in stats)
    has_atac = any(has_atac_layout(st) for st in stats)
    if has_gex and has_atac and not args.prefer_gex:
        raise TenxInferenceError(
            "Both GEX and ATAC barcode evidence was detected. Refusing to mix modalities; rerun with --prefer-gex "
            "to emit only the GEX library for STARsolo, or process ATAC separately with --allow-atac-output."
        )
    if has_gex:
        return classify_gex(stats, args, allow_excluded=args.prefer_gex)
    if has_atac:
        plan = classify_atac(stats, args)
        if args.require_starsolo_compatible or not args.allow_atac_output:
            raise TenxInferenceError(
                "Detected ATAC-only 10x library; it is not a STARsolo Gene Expression input. "
                "Pass --allow-atac-output only when you intentionally want Cell Ranger/ARC ATAC-style FASTQs."
            )
        return plan
    raise TenxInferenceError("No file matched a known 10x whitelist with expected read geometry")

def write_existing_fastq_as_gz(src: Path, dest: Path, action: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if compression_of(src) == "gz":
        if action == "hardlink":
            try:
                os.link(os.path.realpath(src), dest)
                return
            except OSError:
                # Cross-device links are common in Nextflow work dirs; copy is
                # safer than silently creating a symlink when hardlink was asked.
                shutil.copy2(os.path.realpath(src), dest)
                return
        if action == "symlink":
            os.symlink(os.path.realpath(src), dest)
            return
        shutil.copy2(os.path.realpath(src), dest)
        return
    with open_bytes_auto(src) as in_handle, gzip.open(dest, "wb", compresslevel=6) as out_handle:
        shutil.copyfileobj(in_handle, out_handle, length=1024 * 1024)


def write_synthetic_cb_umi(cb_path: Path, umi_path: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(dest, "wt", encoding="utf-8") as out_handle:
        for cb_rec, umi_rec in zip(iter_fastq_records(cb_path), iter_fastq_records(umi_path)):
            cb_h, cb_seq, cb_plus, cb_qual = cb_rec
            umi_h, umi_seq, _umi_plus, umi_qual = umi_rec
            if normalized_read_id(cb_h) != normalized_read_id(umi_h):
                raise TenxInferenceError(
                    f"3prime v1 CB/UMI record ID mismatch at record {n+1}: "
                    f"{normalized_read_id(cb_h)} vs {normalized_read_id(umi_h)}"
                )
            out_handle.write(f"{cb_h}\n{cb_seq}{umi_seq}\n{cb_plus}\n{cb_qual}{umi_qual}\n")
            n += 1
    # Ensure both inputs had the same number of records.
    cb_count = sum(1 for _ in iter_fastq_records(cb_path))
    umi_count = sum(1 for _ in iter_fastq_records(umi_path))
    if cb_count != umi_count or cb_count != n:
        raise TenxInferenceError(f"3prime v1 CB/UMI record counts differ: {cb_count} vs {umi_count}")


def count_fastq_records(path: Path) -> int:
    return sum(1 for _ in iter_fastq_records(path))


def first_normalized_ids(path: Path, n: int) -> List[str]:
    ids: List[str] = []
    for h, _s, _p, _q in iter_fastq_records(path):
        ids.append(normalized_read_id(h))
        if len(ids) >= n:
            break
    return ids


def validate_plan_mates(plan: InferencePlan, args: argparse.Namespace) -> Dict[str, object]:
    checks: Dict[str, object] = {"counts_checked": not args.no_check_counts, "ids_checked": not args.no_check_ids}
    entries_by_lane: Dict[str, List[OutputEntry]] = collections.defaultdict(list)
    for entry in plan.entries:
        entries_by_lane[entry.lane].append(entry)
    if not args.no_check_counts:
        count_report: Dict[str, Dict[str, int]] = {}
        for lane, entries in sorted(entries_by_lane.items()):
            counts: Dict[str, int] = {}
            for entry in entries:
                for src in entry.sources:
                    counts[Path(src).name] = count_fastq_records(Path(src))
            if len(set(counts.values())) > 1:
                raise TenxInferenceError(f"Read counts differ across mates for {lane}: {counts}")
            count_report[lane] = counts
        checks["read_counts"] = count_report
    if not args.no_check_ids:
        id_report: Dict[str, object] = {}
        for lane, entries in sorted(entries_by_lane.items()):
            r1_entries = [e for e in entries if e.read_type == "R1"]
            if not r1_entries:
                continue
            ref_src = Path(r1_entries[0].sources[0])
            ref_ids = first_normalized_ids(ref_src, args.id_check_records)
            lane_report: Dict[str, int] = {}
            for entry in entries:
                for src in entry.sources:
                    p = Path(src)
                    if p == ref_src:
                        continue
                    ids = first_normalized_ids(p, args.id_check_records)
                    n = min(len(ref_ids), len(ids))
                    mismatches = sum(1 for a, b in zip(ref_ids[:n], ids[:n]) if a != b)
                    lane_report[p.name] = mismatches
                    if mismatches:
                        raise TenxInferenceError(
                            f"FASTQ record IDs are not concordant for {lane}: "
                            f"{mismatches}/{n} mismatch between {ref_src.name} and {p.name}"
                        )
            id_report[lane] = lane_report
        checks["id_mismatches"] = id_report
    return checks


def execute_plan(plan: InferencePlan, outdir: Path, action: str, dry_run: bool) -> None:
    if dry_run:
        return
    outdir.mkdir(parents=True, exist_ok=True)
    names = [entry.output_name for entry in plan.entries]
    if len(names) != len(set(names)):
        raise TenxInferenceError(f"Duplicate output names in plan: {names}")
    for entry in plan.entries:
        dest = outdir / entry.output_name
        if entry.synthetic:
            if len(entry.sources) != 2:
                raise TenxInferenceError(f"Synthetic entry requires CB and UMI sources: {entry}")
            write_synthetic_cb_umi(Path(entry.sources[0]), Path(entry.sources[1]), dest)
        else:
            if len(entry.sources) != 1:
                raise TenxInferenceError(f"Non-synthetic entry requires one source: {entry}")
            write_existing_fastq_as_gz(Path(entry.sources[0]), dest, action)


def write_manifest(plan: InferencePlan, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("sample\tchemistry_id\tmodality\tlane\tread_type\toutput_name\tsynthetic\tsources\n")
        for entry in plan.entries:
            handle.write(
                f"{plan.sample}\t{plan.chemistry_id}\t{plan.modality}\t{entry.lane}\t{entry.read_type}\t"
                f"{entry.output_name}\t{entry.synthetic}\t{';'.join(entry.sources)}\n"
            )


def report_dict(plan: Optional[InferencePlan], error: Optional[str] = None, stats: Optional[List[FastqStats]] = None) -> Dict[str, object]:
    if plan is None:
        return {"ok": False, "error": error, "stats": [st.to_json() for st in (stats or [])]}
    r1_umi = next((e.umi_len for e in plan.entries if e.read_type == "R1" and e.umi_len is not None), None)
    cb_len = None
    for st in plan.stats:
        m = best_match(st, plan.modality if plan.modality in {"gex", "atac"} else None)
        if m is not None and ((plan.modality == "gex" and has_gex_layout(st)) or (plan.modality == "atac" and has_atac_layout(st))):
            cb_len = m.cb_len
            break
    return {
        "ok": True,
        "sample": plan.sample,
        "chemistry_id": plan.chemistry_id,
        "chemistry_label": plan.chemistry_label,
        "modality": plan.modality,
        "cb_len": cb_len,
        "umi_len": r1_umi,
        "usable_for_starsolo": plan.usable_for_starsolo,
        "confidence": plan.confidence,
        "warnings": plan.warnings,
        "excluded_files": plan.excluded_files,
        "checks": plan.checks,
        "outputs": [asdict(e) for e in plan.entries],
        "stats": [st.to_json() for st in plan.stats],
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastqs", nargs="*", type=Path, help="Input FASTQ/FASTQ.GZ/FASTQ.BZ2 files for one sample")
    parser.add_argument("--fastqs", dest="fastqs_opt", nargs="+", type=Path, help="Input FASTQ files; compatibility alias for positional inputs")
    parser.add_argument("--sample", "--sample-id", dest="sample", required=True, help="Output sample name used in Cell Ranger FASTQ names")
    parser.add_argument("--whitelists", "--whitelist-dir", dest="whitelists", required=True, type=Path, help="Directory containing 10x barcode whitelist files")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for normalized FASTQs")
    parser.add_argument("--json", dest="json_path", type=Path, default=Path("tenx_fastq_inference.json"))
    parser.add_argument("--manifest", type=Path, default=Path("tenx_fastq_manifest.tsv"))
    parser.add_argument("--lane", type=int, default=None, help="Force a single output lane number, e.g. 1 -> L001")
    parser.add_argument("--sample-records", "--n-sample", dest="sample_records", type=int, default=200000, help="Records sampled per FASTQ")
    parser.add_argument("--min-whitelist-hits", type=int, default=1000)
    parser.add_argument("--min-whitelist-fraction", "--min-frac", dest="min_whitelist_fraction", type=float, default=0.20)
    parser.add_argument("--min-bio-read-length", "--min-cdna", dest="min_bio_read_length", type=int, default=40)
    parser.add_argument("--min-gex-r2-length", type=int, default=60)
    parser.add_argument("--min-atac-genomic-read-length", type=int, default=40, help="Minimum length for ATAC genomic R1/R2 reads")
    parser.add_argument("--allow-v1-umi-lengths", default="5,10", help="Comma-separated 3prime v1 UMI lengths to accept")
    parser.add_argument("--prefer-gex", action="store_true", help="If GEX and ATAC are both detected, emit only GEX inputs")
    parser.add_argument("--allow-atac-output", action="store_true", help="Permit ATAC-only Cell Ranger/ARC-style output; otherwise ATAC-only inputs fail safely")
    parser.add_argument("--require-starsolo-compatible", "--require-gex", action="store_true", help="Fail if the detected library is not GEX/STARsolo compatible")
    parser.add_argument("--allow-missing-whitelists", action="store_true", help="Do not require every known whitelist family to be present; useful only for tests")
    parser.add_argument("--no-check-counts", action="store_true", help="Skip exact read-count equality checks across emitted mates")
    parser.add_argument("--no-check-ids", action="store_true", help="Skip FASTQ record ID concordance checks across emitted mates")
    parser.add_argument("--id-check-records", type=int, default=1000, help="Number of leading records for mate ID concordance checks")
    parser.add_argument(
        "--action",
        "--mode",
        dest="action",
        choices=("copy", "link", "symlink", "hardlink"),
        default="symlink",
        help="How to handle existing gzipped inputs",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.fastqs_opt:
        if args.fastqs:
            parser.error("Use either positional FASTQs or --fastqs, not both")
        args.fastqs = args.fastqs_opt
    delattr(args, "fastqs_opt")
    if not args.fastqs:
        parser.error("At least two FASTQ files are required")
    if len(args.fastqs) < 2:
        parser.error("Need at least two FASTQ files: a barcode read and a biological read")
    if not SAMPLE_RE.fullmatch(args.sample):
        parser.error(f"Unsafe sample name for Cell Ranger FASTQ naming: {args.sample!r}")
    if args.lane is not None and not (1 <= args.lane <= 999):
        parser.error("--lane must be in the range 1..999")
    if args.action == "link":
        args.action = "symlink"
    args.allow_v1_umi_lengths = {int(x) for x in str(args.allow_v1_umi_lengths).split(",") if x.strip()}
    for path in args.fastqs:
        if not path.exists():
            parser.error(f"Input FASTQ does not exist: {path}")
    if not args.whitelists.exists():
        parser.error(f"Whitelist directory does not exist: {args.whitelists}")
    return args

def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    stats: List[FastqStats] = []
    plan: Optional[InferencePlan] = None
    try:
        stats = collect_stats(args.fastqs, args.sample_records)
        match_whitelists(stats, WhitelistStore(args.whitelists), args)
        plan = infer_plan(stats, args)
        apply_forced_lane(plan, args.lane)
        plan.checks = validate_plan_mates(plan, args)
        execute_plan(plan, args.outdir, args.action, args.dry_run)
        write_manifest(plan, args.manifest)
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json_path, "w", encoding="utf-8") as handle:
            json.dump(report_dict(plan), handle, indent=2, sort_keys=True)
        if plan.warnings:
            for warning in plan.warnings:
                print(f"WARNING: {warning}", file=sys.stderr)
        if plan.excluded_files:
            print("Excluded non-selected files: " + ", ".join(plan.excluded_files), file=sys.stderr)
        print(
            f"Detected {plan.chemistry_label} ({plan.chemistry_id}); modality={plan.modality}; "
            f"outputs={len(plan.entries)}",
            file=sys.stderr,
        )
        return 0
    except Exception as exc:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json_path, "w", encoding="utf-8") as handle:
            json.dump(report_dict(None, error=str(exc), stats=stats), handle, indent=2, sort_keys=True)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
