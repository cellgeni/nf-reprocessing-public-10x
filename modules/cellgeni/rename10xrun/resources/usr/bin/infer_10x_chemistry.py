#!/usr/bin/env python3
# =============================================================================
# infer_10x_chemistry.py
# -----------------------------------------------------------------------------
# Infer the 10x Genomics chemistry of a set of FASTQ files by matching reads
# against barcode whitelists, identify the role of each read file (R1 barcode,
# R2 cDNA, I1/I2 sample index, and -- for 3' v1 -- the separate UMI read),
# run a battery of validation checks, and rename the files to the Cell Ranger
# naming convention so they can be fed to Cell Ranger or STARsolo.
#
#   [Sample]_S1_L00[Lane]_[R1|R2|I1|I2]_001.fastq.gz
#
# Design goal: be correct rather than clever. The barcode read is found by
# *whitelist content*, never by length alone, so the tool is robust to
# over-sequenced R1, swapped mates, and mislabelled files. See the companion
# docs/10x_read_geometry.md for the full geometry reference and check list.
#
# This is a single-mate-set tool: one invocation handles one sample / one lane
# (i.e. one SRR / one GSM / one SRS-ERS, matching the reprocessing pipeline's
# "one 10x run per sample" assumption). For multi-lane data, run once per lane.
# =============================================================================

from __future__ import annotations

import argparse
import bz2
import gzip
import json
import os
import random
import re
import shutil
import sys
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Chemistry reference table
# ---------------------------------------------------------------------------
# cb_prefix_len is how many bp of the barcode read we test against the
# whitelist (== cell-barcode length). strand is the STARsolo --soloStrand
# default; "Unknown" means it must be resolved downstream (the 3'v2 / 5'
# ambiguity). split_cb_umi is True only for the original 3' v1 layout.

@dataclass(frozen=True)
class Chemistry:
    key: str
    label: str
    whitelist: str
    cb_len: int
    umi_len: int
    cb_prefix_len: int
    strand: str
    generation: str
    split_cb_umi: bool = False
    notes: str = ""


# 3' v1 was shipped with two UMI lengths over its lifetime: the common 10 bp
# UMI and an early 5 bp UMI used by some datasets (e.g. SRR10759480). The UMI
# read is therefore MEASURED at run time, not assumed -- nominal umi_len below
# is just the common case. See cellgeni/reprocess_public_10x#17.
V1_UMI_LENS = (5, 10)

CHEMISTRIES: list[Chemistry] = [
    Chemistry("3pv1", "3' v1", "737K-april-2014_rc.txt", 14, 10, 14,
              "Forward", "GemCode / Chromium v1", split_cb_umi=True,
              notes="CB (14bp) and UMI (5 or 10bp) sequenced as SEPARATE "
                    "reads; R1 is reconstructed as CB+UMI and the UMI length "
                    "is measured from the data."),
    Chemistry("3pv2_or_5p", "3' v2  or  5' v1/v2", "737K-august-2016.txt", 16, 10, 16,
              "Unknown", "v2 / 5' v1-v2",
              notes="Whitelist shared by 3' v2 and 5' v1/v2. Renaming is "
                    "identical for both; resolve 3' vs 5' (strand) downstream."),
    Chemistry("3pv3", "3' v3 / v3.1", "3M-february-2018.txt", 16, 12, 16,
              "Forward", "v3 / v3.1"),
    Chemistry("3pv4", "3' v4 (GEM-X)", "3M-3pgex-may-2023.txt", 16, 12, 16,
              "Forward", "GEM-X 3'"),
    Chemistry("5pv3", "5' v3 (GEM-X)", "3M-5pgex-jan-2023.txt", 16, 12, 16,
              "Reverse", "GEM-X 5'"),
    Chemistry("multiome", "Multiome GEX (ARC v1)", "737K-arc-v1.txt", 16, 12, 16,
              "Forward", "Multiome ARC v1"),
]
CHEM_BY_KEY = {c.key: c for c in CHEMISTRIES}


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log(level: str, msg: str) -> None:
    sys.stderr.write(f"[{level:<5}] {msg}\n")


def info(msg): log("INFO", msg)
def warn(msg): log("WARN", msg)
def err(msg):  log("ERROR", msg)


class DetectionError(Exception):
    """Raised on a hard validation failure (the tool should exit non-zero)."""


# ---------------------------------------------------------------------------
# Compression helpers (detect by magic bytes, not by extension)
# ---------------------------------------------------------------------------

def compression_of(path: str) -> str:
    with open(path, "rb") as fh:
        magic = fh.read(3)
    if magic[:2] == b"\x1f\x8b":
        return "gz"
    if magic[:3] == b"BZh":
        return "bz2"
    return "none"


def open_text(path: str):
    comp = compression_of(path)
    if comp == "gz":
        return gzip.open(path, "rt")
    if comp == "bz2":
        return bz2.open(path, "rt")
    return open(path, "rt")


def open_binary(path: str):
    comp = compression_of(path)
    if comp == "gz":
        return gzip.open(path, "rb")
    if comp == "bz2":
        return bz2.open(path, "rb")
    return open(path, "rb")


# ---------------------------------------------------------------------------
# FASTQ sampling and read-level stats
# ---------------------------------------------------------------------------

@dataclass
class FastqInfo:
    path: str
    compression: str
    n_seen: int                 # reads scanned
    median_len: int             # typical read length (from the sample)
    min_len: int
    max_len: int
    distinct_len: int           # number of distinct lengths seen in the sample
    prefix14: list[str] = field(default_factory=list)
    prefix16: list[str] = field(default_factory=list)
    # filled in later:
    role: Optional[str] = None
    wl_match: dict[str, float] = field(default_factory=dict)


def sample_fastq(path: str, n_sample: int, n_scan: int, seed: int = 100) -> FastqInfo:
    """Reservoir-sample up to n_sample reads from the first n_scan reads.

    Stats (length median/min/max, distinct lengths) are computed from the
    reservoir, which is a uniform sample of the scanned reads. We cap the scan
    at n_scan to bound runtime on very large files; the first ~2M reads are
    plenty representative for chemistry detection.
    """
    rng = random.Random(seed)
    res: list[str] = []
    n_seen = 0
    with open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline()
            fh.readline()          # plus line
            qual = fh.readline()
            if not qual:
                break              # truncated final record
            seq = seq.rstrip("\n")
            n_seen += 1
            if len(res) < n_sample:
                res.append(seq)
            else:
                j = rng.randint(0, n_seen - 1)
                if j < n_sample:
                    res[j] = seq
            if n_scan and n_seen >= n_scan:
                break

    if not res:
        raise DetectionError(f"No reads found in {path}")

    lengths = sorted(len(s) for s in res)
    median_len = lengths[len(lengths) // 2]
    return FastqInfo(
        path=path,
        compression=compression_of(path),
        n_seen=n_seen,
        median_len=median_len,
        min_len=lengths[0],
        max_len=lengths[-1],
        distinct_len=len(set(lengths)),
        prefix14=[s[:14] for s in res],
        prefix16=[s[:16] for s in res],
    )


def count_reads(path: str) -> int:
    """Exact read count by counting newlines (4 lines per record)."""
    n = 0
    with open_binary(path) as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            n += chunk.count(b"\n")
    return n // 4


# ---------------------------------------------------------------------------
# Whitelist loading and matching
# ---------------------------------------------------------------------------

def load_whitelist(path: str, k: int) -> set[str]:
    """Load the first-k-bp prefix of each whitelist barcode into a set.

    Robust to translation-style whitelists (e.g. 3M-5pgex-jan-2023): we take
    the first whitespace-delimited field and its first k characters, which is
    the cell barcode for all Cell Ranger whitelist files.
    """
    s: set[str] = set()
    with open_text(path) as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            bc = tok[0]
            if len(bc) >= k:
                s.add(bc[:k])
    return s


def match_against(infos: list[FastqInfo], wl: set[str], k: int) -> dict[str, float]:
    out: dict[str, float] = {}
    for fi in infos:
        pref = fi.prefix14 if k == 14 else fi.prefix16
        if not pref:
            out[fi.path] = 0.0
            continue
        # only count prefixes that are full length k (short reads -> not a match)
        hits = sum(1 for p in pref if len(p) == k and p in wl)
        out[fi.path] = hits / len(pref)
    return out


def detect_chemistry(infos: list[FastqInfo], wl_dir: str, min_frac: float):
    """Return (chemistry, barcode_path). Populates fi.wl_match for every file.

    Whitelists are loaded one at a time to keep peak memory ~ one whitelist
    (the 3M files are ~3M lines).
    """
    present = []
    for c in CHEMISTRIES:
        wlpath = os.path.join(wl_dir, c.whitelist)
        if os.path.exists(wlpath):
            present.append((c, wlpath))
        else:
            warn(f"whitelist not found, skipping: {wlpath}")
    if not present:
        raise DetectionError(f"No whitelists found in {wl_dir}")

    for c, wlpath in present:
        wl = load_whitelist(wlpath, c.cb_prefix_len)
        info(f"  loaded {c.whitelist} ({len(wl):,} barcodes)")
        frac = match_against(infos, wl, c.cb_prefix_len)
        for fi in infos:
            fi.wl_match[c.key] = frac[fi.path]
        del wl

    # Collect strong matches across all (file, chemistry) pairs.
    strong = [
        (fi, ckey, f)
        for fi in infos
        for ckey, f in fi.wl_match.items()
        if f >= min_frac
    ]
    if not strong:
        best = max(
            ((fi, ck, f) for fi in infos for ck, f in fi.wl_match.items()),
            key=lambda t: t[2],
        )
        raise DetectionError(
            "No read file matched any 10x whitelist above the threshold "
            f"({min_frac:.0%}). Best was {os.path.basename(best[0].path)} "
            f"vs {CHEM_BY_KEY[best[1]].label} at {best[2]:.1%}. "
            "This is most likely not 10x single-cell data."
        )

    distinct_chems = {ckey for _, ckey, _ in strong}
    if len(distinct_chems) > 1:
        detail = ", ".join(
            f"{os.path.basename(fi.path)}~{CHEM_BY_KEY[ck].label}={f:.0%}"
            for fi, ck, f in strong
        )
        raise DetectionError(
            f"Ambiguous chemistry: more than one whitelist matched strongly "
            f"({detail}). Refusing to guess."
        )

    ckey = distinct_chems.pop()
    chem = CHEM_BY_KEY[ckey]

    barcode_hits = [(fi, f) for fi, ck, f in strong if ck == ckey]
    if len(barcode_hits) > 1:
        detail = ", ".join(f"{os.path.basename(fi.path)}={f:.0%}" for fi, f in barcode_hits)
        raise DetectionError(
            f"More than one file matched the {chem.label} whitelist "
            f"({detail}). Cannot identify a unique barcode read."
        )

    barcode_fi = barcode_hits[0][0]
    info(f"Detected chemistry: {chem.label}  "
         f"(barcode read = {os.path.basename(barcode_fi.path)}, "
         f"match = {barcode_hits[0][1]:.1%})")
    return chem, barcode_fi


# ---------------------------------------------------------------------------
# Role assignment
# ---------------------------------------------------------------------------

_INDEX_HINT = re.compile(r"_(I[12])[._]")
_READ_HINT = re.compile(r"_(R[123])[._]")


def filename_index_hint(path: str) -> Optional[str]:
    m = _INDEX_HINT.search(os.path.basename(path))
    return m.group(1) if m else None


@dataclass
class Roles:
    chem: Chemistry
    barcode: FastqInfo                       # -> R1 (or CB part for v1)
    cdna: FastqInfo                          # -> R2
    umi: Optional[FastqInfo] = None          # v1 only -> merged into R1
    index1: Optional[FastqInfo] = None       # -> I1
    index2: Optional[FastqInfo] = None       # -> I2
    index_order_confident: bool = True
    effective_umi_len: int = 0
    warnings: list[str] = field(default_factory=list)


def assign_roles(infos: list[FastqInfo], chem: Chemistry, barcode_fi: FastqInfo,
                 min_cdna: int, max_index: int) -> Roles:
    others = [fi for fi in infos if fi is not barcode_fi]
    roles = Roles(chem=chem, barcode=barcode_fi, cdna=None)  # type: ignore
    barcode_fi.role = "R1"

    umi_fi = None
    if chem.split_cb_umi:
        # 3' v1: the cell barcode (14 bp) and the UMI are sequenced as TWO
        # separate reads, and the UMI is either 5 or 10 bp depending on the
        # run -- some early v1 datasets use a 5 bp UMI (e.g. SRR10759480; see
        # cellgeni/reprocess_public_10x#17). We MEASURE it, never assume 10.
        #
        # Set the long cDNA read aside first; the UMI is then the short read
        # that remains. If a sample-index read is also present we may see two
        # short reads, so choose the one nearest a valid v1 UMI length and,
        # on a tie, the one that does NOT carry an I1/I2 filename hint.
        shorts = [fi for fi in others if fi.median_len < min_cdna]
        if not shorts:
            raise DetectionError(
                f"Detected a 3' v1 cell-barcode read but found no separate "
                f"short UMI read (< {min_cdna} bp). The SRA dump probably did "
                "not preserve the UMI read, so CB+UMI cannot be reconstructed."
            )

        def _umi_rank(fi: FastqInfo):
            near = min(abs(fi.median_len - L) for L in V1_UMI_LENS)
            has_index_hint = filename_index_hint(fi.path) is not None
            return (near, has_index_hint, fi.median_len)

        shorts.sort(key=_umi_rank)
        umi_fi = shorts[0]
        if len(shorts) > 1:
            roles.warnings.append(
                "More than one short read accompanies the v1 barcode read; "
                f"chose {os.path.basename(umi_fi.path)} ({umi_fi.median_len} bp) "
                "as the UMI by closeness to a valid v1 UMI length (5 or 10 bp). "
                "If a sample-index read is also present, double-check this pick."
            )
        if min(abs(umi_fi.median_len - L) for L in V1_UMI_LENS) > 2:
            roles.warnings.append(
                f"v1 UMI read is {umi_fi.median_len} bp, which is not close to "
                "either documented v1 UMI length (5 or 10 bp); proceeding with "
                "the measured length, but please verify the file roles."
            )
        umi_fi.role = "UMI"
        roles.umi = umi_fi
        others = [fi for fi in others if fi is not umi_fi]
        roles.effective_umi_len = umi_fi.median_len
    else:
        roles.effective_umi_len = chem.umi_len

    # cDNA = longest remaining read that clears the cDNA length floor.
    cdna_cands = sorted(others, key=lambda fi: fi.median_len, reverse=True)
    cdna_fi = next((fi for fi in cdna_cands if fi.median_len >= min_cdna), None)
    if cdna_fi is None:
        longest = cdna_cands[0] if cdna_cands else None
        extra = f" (longest was {longest.median_len} bp)" if longest else ""
        raise DetectionError(
            f"No cDNA (R2) read >= {min_cdna} bp found{extra}. "
            "If this is multiome ATAC data it is out of scope for this module."
        )
    cdna_fi.role = "R2"
    roles.cdna = cdna_fi
    others = [fi for fi in others if fi is not cdna_fi]

    # Remaining files: sample index reads (I1/I2), expected to be short.
    index_fis = []
    for fi in others:
        if fi.median_len <= max_index:
            index_fis.append(fi)
        else:
            roles.warnings.append(
                f"Unexpected extra read {os.path.basename(fi.path)} "
                f"({fi.median_len} bp) is neither cDNA nor a short index; ignoring."
            )
    if len(index_fis) > 2:
        roles.warnings.append(
            f"Found {len(index_fis)} candidate index reads; using the two shortest."
        )
        index_fis = sorted(index_fis, key=lambda fi: fi.median_len)[:2]

    # Order I1/I2: trust an explicit filename hint, else fall back to encounter
    # order and flag that we are not sure which is i7 vs i5.
    if index_fis:
        hinted = {filename_index_hint(fi.path): fi for fi in index_fis}
        hinted.pop(None, None)
        if "I1" in hinted or "I2" in hinted:
            roles.index1 = hinted.get("I1")
            roles.index2 = hinted.get("I2")
            leftovers = [fi for fi in index_fis if fi not in (roles.index1, roles.index2)]
            if roles.index1 is None and leftovers:
                roles.index1 = leftovers.pop(0)
            if roles.index2 is None and leftovers:
                roles.index2 = leftovers.pop(0)
        else:
            roles.index1 = index_fis[0]
            roles.index2 = index_fis[1] if len(index_fis) > 1 else None
            if roles.index2 is not None:
                roles.index_order_confident = False
                roles.warnings.append(
                    "Two index reads present but neither carries an I1/I2 "
                    "filename hint; I1/I2 assignment is arbitrary. This does "
                    "not affect STARsolo or Cell Ranger reprocessing."
                )
        for fi, lab in ((roles.index1, "I1"), (roles.index2, "I2")):
            if fi is not None:
                fi.role = lab

    return roles


# ---------------------------------------------------------------------------
# Validation checks
# ---------------------------------------------------------------------------

def normalize_id(header: str) -> str:
    x = header[1:] if header.startswith("@") else header
    x = x.split()[0] if x.split() else x
    return re.sub(r"[./][123]$", "", x)


def first_ids(path: str, n: int) -> list[str]:
    out = []
    with open_text(path) as fh:
        while len(out) < n:
            h = fh.readline()
            if not h:
                break
            fh.readline(); fh.readline(); fh.readline()
            out.append(normalize_id(h.rstrip("\n")))
    return out


def validate(roles: Roles, check_counts: bool, check_ids: bool) -> dict:
    chem = roles.chem
    report: dict[str, object] = {}

    # C2: barcode read long enough for CB (+UMI). Mirror cellgeni: shrink UMI
    # if R1 is a little short; fatal only if it cannot hold the cell barcode.
    if chem.split_cb_umi:
        if roles.barcode.median_len < chem.cb_len:
            raise DetectionError(
                f"v1 CB read is {roles.barcode.median_len} bp, "
                f"shorter than the {chem.cb_len} bp cell barcode.")
    else:
        bc_umi = chem.cb_len + chem.umi_len
        r1 = roles.barcode.median_len
        if r1 < chem.cb_len:
            raise DetectionError(
                f"Barcode read (R1) is {r1} bp, shorter than the "
                f"{chem.cb_len} bp cell barcode.")
        if r1 < bc_umi:
            new_umi = r1 - chem.cb_len
            warn(f"R1 ({r1} bp) < CB+UMI ({bc_umi} bp); "
                 f"effective UMI length reduced to {new_umi}.")
            roles.effective_umi_len = new_umi
    report["C2_barcode_length_ok"] = True

    # C3: barcode read should be (near) fixed length. A short, ragged R1 is the
    # classic sign of a quality-trimmed barcode and is not safe to use.
    if roles.barcode.distinct_len > 1 and roles.barcode.median_len <= 30:
        raise DetectionError(
            f"Barcode read has {roles.barcode.distinct_len} distinct lengths "
            f"around {roles.barcode.median_len} bp; looks quality-trimmed.")
    report["C3_barcode_fixed_length"] = roles.barcode.distinct_len == 1

    # C4: cDNA length sanity (warn-level below 40, already fatal if absent).
    if roles.cdna.median_len < 40:
        warn(f"cDNA read (R2) is only {roles.cdna.median_len} bp; "
             "transcriptome alignment rate may suffer.")
    report["C4_cdna_length"] = roles.cdna.median_len

    # C5: equal read counts across every mate that will be emitted.
    if check_counts:
        members = {"R1_src": roles.barcode, "R2": roles.cdna}
        if roles.umi:
            members["UMI"] = roles.umi
        if roles.index1:
            members["I1"] = roles.index1
        if roles.index2:
            members["I2"] = roles.index2
        counts = {name: count_reads(fi.path) for name, fi in members.items()}
        report["C5_read_counts"] = counts
        uniq = set(counts.values())
        if len(uniq) > 1:
            raise DetectionError(f"Read counts differ across mates: {counts}")
        info(f"Read count per mate: {next(iter(uniq)):,} (consistent)")
    else:
        report["C5_read_counts"] = "skipped"

    # C8: light record-ID concordance check on the first reads.
    if check_ids:
        ref = first_ids(roles.barcode.path, 1000)
        for name, fi in (("R2", roles.cdna),
                         ("UMI", roles.umi),
                         ("I1", roles.index1),
                         ("I2", roles.index2)):
            if fi is None:
                continue
            other = first_ids(fi.path, 1000)
            n = min(len(ref), len(other))
            mism = sum(1 for a, b in zip(ref[:n], other[:n]) if a != b)
            if mism:
                roles.warnings.append(
                    f"{mism}/{n} record IDs differ between R1 and {name} "
                    "in the first reads (could be a harmless ID-format "
                    "difference, or the files may not be mates).")
        report["C8_id_concordance_checked"] = True

    return report


# ---------------------------------------------------------------------------
# Emitting renamed / reconstructed files
# ---------------------------------------------------------------------------

def cr_name(sample: str, lane: int, read_type: str) -> str:
    return f"{sample}_S1_L{lane:03d}_{read_type}_001.fastq.gz"


def place_gz(src: str, dst: str, mode: str) -> None:
    """Put src at dst as a .gz file.

    gz source: hardlink / symlink / copy according to mode.
    bz2 or plain source: stream-recompress to gz (Cell Ranger cannot read bz2,
    and the .gz name must not lie about the contents).
    """
    comp = compression_of(src)
    if comp == "gz":
        if os.path.lexists(dst):
            os.remove(dst)
        if mode == "hardlink":
            try:
                os.link(os.path.realpath(src), dst)
                return
            except OSError:
                pass  # cross-device etc. -> fall through
        if mode in ("hardlink", "symlink"):
            os.symlink(os.path.realpath(src), dst)
            return
        shutil.copy2(os.path.realpath(src), dst)
        return
    # recompress
    info(f"  recompressing {os.path.basename(src)} ({comp}) -> gzip")
    with open_binary(src) as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)


def merge_cb_umi(cb_path: str, umi_path: str, dst: str, cb_len: int) -> None:
    """Reconstruct a STARsolo-compatible R1 = CB + UMI for 3' v1 data."""
    info(f"  reconstructing R1 = CB({cb_len}) + UMI from "
         f"{os.path.basename(cb_path)} + {os.path.basename(umi_path)}")
    n = 0
    with open_text(cb_path) as fcb, open_text(umi_path) as fumi, \
            gzip.open(dst, "wt") as out:
        while True:
            h1 = fcb.readline()
            if not h1:
                break
            s1 = fcb.readline().rstrip("\n"); fcb.readline(); q1 = fcb.readline().rstrip("\n")
            h2 = fumi.readline()
            if not h2:
                raise DetectionError("UMI read ended before barcode read during merge.")
            s2 = fumi.readline().rstrip("\n"); fumi.readline(); q2 = fumi.readline().rstrip("\n")
            if normalize_id(h1.rstrip("\n")) != normalize_id(h2.rstrip("\n")):
                raise DetectionError(
                    f"CB/UMI record mismatch at read {n+1}: "
                    f"{h1.strip()} vs {h2.strip()}")
            out.write(f"{h1.rstrip(chr(10))}\n{s1[:cb_len]}{s2}\n+\n{q1[:cb_len]}{q2}\n")
            n += 1
    info(f"  wrote {n:,} reconstructed R1 records")


def emit(roles: Roles, sample: str, lane: int, outdir: str, mode: str) -> dict:
    os.makedirs(outdir, exist_ok=True)
    written: dict[str, str] = {}

    if roles.chem.split_cb_umi:
        r1 = os.path.join(outdir, cr_name(sample, lane, "R1"))
        merge_cb_umi(roles.barcode.path, roles.umi.path, r1, cb_len=roles.chem.cb_len)
        written["R1"] = r1
    else:
        r1 = os.path.join(outdir, cr_name(sample, lane, "R1"))
        place_gz(roles.barcode.path, r1, mode)
        written["R1"] = r1

    r2 = os.path.join(outdir, cr_name(sample, lane, "R2"))
    place_gz(roles.cdna.path, r2, mode)
    written["R2"] = r2

    for fi, lab in ((roles.index1, "I1"), (roles.index2, "I2")):
        if fi is None:
            continue
        dst = os.path.join(outdir, cr_name(sample, lane, lab))
        place_gz(fi.path, dst, mode)
        written[lab] = dst

    for lab, p in written.items():
        info(f"  {lab} -> {os.path.basename(p)}")
    return written


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_manifest(sample, lane, roles, written, checks, args) -> dict:
    chem = roles.chem
    def src(fi): return os.path.basename(fi.path) if fi else None
    return {
        "sample": sample,
        "lane": lane,
        "chemistry": chem.key,
        "chemistry_label": chem.label,
        "generation": chem.generation,
        "whitelist": chem.whitelist,
        "cb_length": chem.cb_len,
        "umi_length": roles.effective_umi_len,
        "umi_length_nominal": chem.umi_len,
        "strand": chem.strand,                       # "Unknown" -> resolve downstream
        "strand_resolved_downstream": chem.strand == "Unknown",
        "split_cb_umi": chem.split_cb_umi,
        "index_order_confident": roles.index_order_confident,
        "file_roles": {
            "R1_source": src(roles.barcode),
            "UMI_source": src(roles.umi),
            "R2_source": src(roles.cdna),
            "I1_source": src(roles.index1),
            "I2_source": src(roles.index2),
        },
        "renamed": {k: os.path.basename(v) for k, v in written.items()},
        "whitelist_match_fraction": {
            os.path.basename(fi.path): {k: round(v, 4) for k, v in fi.wl_match.items()}
            for fi in roles_all_infos(roles)
        },
        "read_geometry": {
            os.path.basename(fi.path): {
                "median_len": fi.median_len,
                "min_len": fi.min_len,
                "max_len": fi.max_len,
                "distinct_len": fi.distinct_len,
                "role": fi.role,
            } for fi in roles_all_infos(roles)
        },
        "checks": checks,
        "warnings": roles.warnings,
        "tool": {
            "name": "infer_10x_chemistry.py",
            "n_sample": args.n_sample,
            "n_scan": args.n_scan,
            "min_frac": args.min_frac,
        },
    }


def roles_all_infos(roles: Roles) -> list[FastqInfo]:
    out = [roles.barcode, roles.cdna]
    for fi in (roles.umi, roles.index1, roles.index2):
        if fi is not None:
            out.append(fi)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Infer 10x chemistry from FASTQ geometry + whitelists and "
                    "rename to the Cell Ranger convention.")
    p.add_argument("--fastqs", nargs="+", required=True,
                   help="All FASTQ files for ONE sample/lane (2-4 files).")
    p.add_argument("--sample-id", required=True, help="Output sample name.")
    p.add_argument("--whitelist-dir",
                   default=os.environ.get("TENX_WHITELIST_DIR",
                                          "/nfs/cellgeni/STAR/whitelists"),
                   help="Directory holding the 10x barcode whitelist .txt files.")
    p.add_argument("--outdir", default="renamed", help="Where to write renamed files.")
    p.add_argument("--lane", type=int, default=1, help="Lane number (default 1).")
    p.add_argument("--json", default=None,
                   help="Path for the chemistry manifest JSON "
                        "(default <outdir>/<sample>.chemistry.json).")
    p.add_argument("--mode", choices=["hardlink", "symlink", "copy"],
                   default="hardlink",
                   help="How to place already-gzipped reads (default hardlink).")
    p.add_argument("--n-sample", type=int, default=200_000,
                   help="Reads to sample for detection (default 200000).")
    p.add_argument("--n-scan", type=int, default=2_000_000,
                   help="Cap on reads scanned per file (default 2000000).")
    p.add_argument("--min-frac", type=float, default=0.25,
                   help="Min fraction of sampled reads matching a whitelist "
                        "to call the barcode read (default 0.25).")
    p.add_argument("--min-cdna", type=int, default=40,
                   help="Minimum median length to accept a read as cDNA/R2.")
    p.add_argument("--max-index", type=int, default=14,
                   help="Maximum median length for a read to count as an index.")
    p.add_argument("--no-check-counts", action="store_true",
                   help="Skip the exact read-count equality check (faster).")
    p.add_argument("--no-check-ids", action="store_true",
                   help="Skip the record-ID concordance check.")
    args = p.parse_args(argv)

    fastqs = []
    for f in args.fastqs:
        if not os.path.exists(f):
            err(f"file not found: {f}")
            return 2
        fastqs.append(f)
    if len(fastqs) < 2:
        err("Need at least 2 FASTQ files (a barcode read and a cDNA read).")
        return 2

    try:
        info(f"Sampling {len(fastqs)} FASTQ files "
             f"(n_sample={args.n_sample:,}, n_scan={args.n_scan:,}) ...")
        infos = [sample_fastq(f, args.n_sample, args.n_scan) for f in fastqs]
        for fi in infos:
            info(f"  {os.path.basename(fi.path):<40} "
                 f"median={fi.median_len:>4}bp  range={fi.min_len}-{fi.max_len}  "
                 f"distinct_len={fi.distinct_len}  comp={fi.compression}")

        info("Matching reads against whitelists ...")
        chem, barcode_fi = detect_chemistry(infos, args.whitelist_dir, args.min_frac)

        roles = assign_roles(infos, chem, barcode_fi, args.min_cdna, args.max_index)
        checks = validate(roles, not args.no_check_counts, not args.no_check_ids)

        written = emit(roles, args.sample_id, args.lane, args.outdir, args.mode)
        manifest = build_manifest(args.sample_id, args.lane, roles, written, checks, args)

        json_path = args.json or os.path.join(args.outdir,
                                              f"{args.sample_id}.chemistry.json")
        with open(json_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        info(f"Wrote manifest: {json_path}")

        for w in roles.warnings:
            warn(w)
        info("DONE.")
        return 0

    except DetectionError as e:
        err(str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())
