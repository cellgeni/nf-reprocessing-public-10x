#!/usr/bin/env python3
"""
Infer one 10x run's library layout from 2-4 FASTQ files, validate the mate set,
and emit Cell Ranger-style FASTQ names for that single run.

The intended pipeline placement is immediately after sra2fastq/bam2fastq/ENA
FASTQ download, before runs are collected per sample. A later naming/grouping
step may compute --fastq-prefix, e.g. SAMPLE_S1_L001; this script only validates
one run and applies that prefix to the inferred read roles.
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
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

DNA_RE = re.compile(r"^[ACGTNacgtn.]+$")
EXPLICIT_ROLE_RE = re.compile(r"(?:^|[_\-.])(R[1234]|I[12])(?:[_\-.]|$)", re.IGNORECASE)
ORDINAL_RE = re.compile(r"(?:^|[_\-.])(?:R)?([1234])(?:\.f(?:ast)?q|\.fastq|\.fq|[_\-.]|$)", re.IGNORECASE)
CR_PREFIX_RE = re.compile(r"^.+_S\d+_L\d{3}$")
SAFE_PREFIX_RE = re.compile(r"^[A-Za-z0-9_.-]+_S\d+_L\d{3}$")


@dataclass(frozen=True)
class ChemistryDef:
    id: str
    label: str
    modality: str          # gex or atac
    whitelist_names: Tuple[str, ...]
    cb_len: int
    umi_len: Optional[int]
    version: str
    chemistry_group: str
    strand_hint: str
    split_v1: bool = False
    barcode_offsets: Tuple[int, ...] = (0,)
    notes: str = ""


CHEMISTRIES: Tuple[ChemistryDef, ...] = (
    ChemistryDef(
        id="gex_3pv1",
        label="Chromium/GemCode Single Cell 3' v1",
        modality="gex",
        whitelist_names=("737K-april-2014_rc.txt", "737K-april-2014_rc.txt.gz"),
        cb_len=14,
        umi_len=None,
        version="3pv1",
        chemistry_group="gex_3pv1",
        strand_hint="Forward",
        split_v1=True,
        notes="3' v1 can be split CB(14)+UMI(5/10); R1 is reconstructed as CB+UMI.",
    ),
    ChemistryDef(
        id="gex_3pv2_or_5pv1v2",
        label="Chromium 3' v2 or 5' v1/v2",
        modality="gex",
        whitelist_names=("737K-august-2016.txt", "737K-august-2016.txt.gz"),
        cb_len=16,
        umi_len=10,
        version="3pv2_or_5pv1v2",
        chemistry_group="gex_737K_august_2016_layout",
        strand_hint="Unknown",
        notes="FASTQ layout cannot separate 3' v2 from 5' v1/v2.",
    ),
    ChemistryDef(
        id="gex_3pv3_family",
        label="Chromium Single Cell 3' v3/v3.1/LT/HT family",
        modality="gex",
        whitelist_names=("3M-february-2018.txt", "3M-february-2018.txt.gz"),
        cb_len=16,
        umi_len=12,
        version="3pv3_family",
        chemistry_group="gex_3pv3_family",
        strand_hint="Forward",
        notes="v3, v3.1, LT and HT share the same FASTQ geometry family.",
    ),
    ChemistryDef(
        id="gex_3pv4_gemx",
        label="GEM-X Universal 3' v4",
        modality="gex",
        whitelist_names=("3M-3pgex-may-2023.txt", "3M-3pgex-may-2023.txt.gz"),
        cb_len=16,
        umi_len=12,
        version="3pv4_gemx",
        chemistry_group="gex_3pv4_gemx",
        strand_hint="Forward",
    ),
    ChemistryDef(
        id="gex_5pv3_gemx",
        label="GEM-X Universal 5' v3",
        modality="gex",
        whitelist_names=("3M-5pgex-jan-2023.txt", "3M-5pgex-jan-2023.txt.gz"),
        cb_len=16,
        umi_len=12,
        version="5pv3_gemx",
        chemistry_group="gex_5pv3_gemx",
        strand_hint="Reverse",
    ),
    ChemistryDef(
        id="gex_multiome_arc_v1",
        label="Chromium Single Cell Multiome Gene Expression ARC v1",
        modality="gex",
        whitelist_names=("gex_737K-arc-v1.txt", "gex_737K-arc-v1.txt.gz", "737K-arc-v1.txt", "737K-arc-v1.txt.gz"),
        cb_len=16,
        umi_len=12,
        version="multiome_arc_v1_gex",
        chemistry_group="multiome_arc_v1_gex",
        strand_hint="Forward",
        notes="GEX R1 contains CB+UMI; ATAC barcode-index reads from the same ARC whitelist are rejected by geometry.",
    ),
    ChemistryDef(
        id="atac_multiome_arc_v1",
        label="Chromium Single Cell Multiome ATAC ARC v1",
        modality="atac",
        whitelist_names=("atac_737K-arc-v1.txt", "atac_737K-arc-v1.txt.gz", "737K-arc-v1.txt", "737K-arc-v1.txt.gz"),
        cb_len=16,
        umi_len=None,
        version="multiome_arc_v1_atac",
        chemistry_group="multiome_arc_v1_atac",
        strand_hint="None",
        barcode_offsets=(0, 8),
        notes="ATAC barcode is I2-like: 16 nt, or 24 nt with an 8-cycle spacer/dark segment.",
    ),
    ChemistryDef(
        id="atac_v1_v2",
        label="Chromium Single Cell ATAC v1/v1.1/v2",
        modality="atac",
        whitelist_names=("737K-cratac-v1.txt", "737K-cratac-v1.txt.gz"),
        cb_len=16,
        umi_len=None,
        version="atac_v1_v2",
        chemistry_group="atac_v1_v2",
        strand_hint="None",
        barcode_offsets=(0, 8),
    ),
)
CHEM_BY_ID = {c.id: c for c in CHEMISTRIES}
V1_UMI_LENGTHS_DEFAULT = {5, 10}
ATAC_BARCODE_LENGTHS = {16, 24}
INDEX_LENGTHS = {6, 7, 8, 9, 10, 14, 16, 24}


class TenxRunError(RuntimeError):
    pass


@dataclass
class Match:
    chemistry_id: str
    modality: str
    label: str
    offset: int
    hits: int
    sampled: int
    fraction: float


@dataclass
class FastqInfo:
    path: str
    basename: str
    compression: str
    sampled_records: int
    length_counts: Dict[int, int]
    common_length: int
    common_fraction: float
    min_length: int
    max_length: int
    explicit_role: Optional[str]
    ordinal: Optional[int]
    matches: List[Match] = field(default_factory=list)
    inferred_role: Optional[str] = None

    def to_json(self) -> Dict[str, object]:
        d = asdict(self)
        d["length_counts"] = {str(k): v for k, v in sorted(self.length_counts.items())}
        return d


@dataclass
class OutputEntry:
    read_type: str
    output_name: str
    sources: List[str]
    source_lengths: List[int]
    synthetic: bool = False
    note: str = ""


@dataclass
class Plan:
    ok: bool
    run_id: str
    sample_id: str
    fastq_prefix: str
    modality: str
    chemistry_id: str
    chemistry_label: str
    chemistry_version: str
    chemistry_group: str
    whitelist: str
    cb_len: Optional[int]
    umi_len: Optional[int]
    strand_hint: str
    layout_confidence: str
    read_lengths: Dict[str, int]
    source_role_map: Dict[str, str]
    outputs: List[OutputEntry]
    excluded_files: List[str]
    warnings: List[str]
    stats: List[FastqInfo]
    checks: Dict[str, object]


# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------


def log(level: str, msg: str) -> None:
    print(f"[{level:<5}] {msg}", file=sys.stderr)


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
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    if comp == "bz2":
        return bz2.open(path, "rt", encoding="utf-8", errors="strict")
    return open(path, "rt", encoding="utf-8", errors="strict")


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
                raise TenxRunError(f"{path}: truncated FASTQ record at record {recno}")
            h = h.rstrip("\r\n")
            s = s.rstrip("\r\n")
            p = p.rstrip("\r\n")
            q = q.rstrip("\r\n")
            if not h.startswith("@"):
                raise TenxRunError(f"{path}: FASTQ header at record {recno} does not start with @")
            if not p.startswith("+"):
                raise TenxRunError(f"{path}: FASTQ plus line at record {recno} does not start with +")
            if len(s) != len(q):
                raise TenxRunError(f"{path}: sequence/quality length mismatch at record {recno}")
            if not DNA_RE.match(s):
                raise TenxRunError(f"{path}: non-DNA sequence characters at record {recno}")
            yield h, s.upper(), p, q


def normalize_read_id(header: str) -> str:
    x = header[1:] if header.startswith("@") else header
    x = x.split()[0]
    return re.sub(r"([/.][1234])$", "", x)


def parse_explicit_role(name: str) -> Optional[str]:
    m = EXPLICIT_ROLE_RE.search(name)
    return m.group(1).upper() if m else None


def parse_ordinal(name: str) -> Optional[int]:
    m = ORDINAL_RE.search(name)
    return int(m.group(1)) if m else None


def natural_key(value: object) -> Tuple[object, ...]:
    return tuple(int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value)))


def collect_fastq_info(paths: Sequence[Path], sample_records: int) -> Tuple[List[FastqInfo], Dict[str, List[str]]]:
    infos: List[FastqInfo] = []
    seq_samples: Dict[str, List[str]] = {}
    resolved: set[Path] = set()
    for path in paths:
        real = path.resolve()
        if real in resolved:
            raise TenxRunError(f"Duplicate input FASTQ after resolving links: {path}")
        resolved.add(real)
        counts: collections.Counter[int] = collections.Counter()
        seqs: List[str] = []
        for _h, seq, _p, _q in iter_fastq_records(path):
            counts[len(seq)] += 1
            seqs.append(seq)
            if len(seqs) >= sample_records:
                break
        if not seqs:
            raise TenxRunError(f"{path}: no FASTQ records sampled")
        common_len, common_n = counts.most_common(1)[0]
        info = FastqInfo(
            path=str(path),
            basename=path.name,
            compression=compression_of(path),
            sampled_records=len(seqs),
            length_counts=dict(counts),
            common_length=common_len,
            common_fraction=common_n / len(seqs),
            min_length=min(counts),
            max_length=max(counts),
            explicit_role=parse_explicit_role(path.name),
            ordinal=parse_ordinal(path.name),
        )
        infos.append(info)
        seq_samples[info.path] = seqs
    return infos, seq_samples


def count_records_strict(path: Path) -> int:
    return sum(1 for _ in iter_fastq_records(path))


def first_ids(path: Path, n: int) -> List[str]:
    out: List[str] = []
    for h, _s, _p, _q in iter_fastq_records(path):
        out.append(normalize_read_id(h))
        if len(out) >= n:
            break
    return out


# ---------------------------------------------------------------------------
# Whitelists and matching
# ---------------------------------------------------------------------------


class WhitelistStore:
    def __init__(self, whitelist_dir: Path):
        self.whitelist_dir = whitelist_dir
        self.paths: Dict[str, Optional[Path]] = {c.id: self._find_path(c) for c in CHEMISTRIES}
        self.cache: Dict[str, set[str]] = {}

    def _find_path(self, chem: ChemistryDef) -> Optional[Path]:
        for name in chem.whitelist_names:
            p = self.whitelist_dir / name
            if p.exists():
                return p
        return None

    def available(self) -> List[ChemistryDef]:
        return [c for c in CHEMISTRIES if self.paths[c.id] is not None]

    def missing(self) -> List[str]:
        return [c.id for c in CHEMISTRIES if self.paths[c.id] is None]

    def path_for(self, chem_id: str) -> str:
        p = self.paths.get(chem_id)
        return p.name if p is not None else ""

    def get(self, chem: ChemistryDef) -> set[str]:
        if chem.id in self.cache:
            return self.cache[chem.id]
        path = self.paths.get(chem.id)
        if path is None:
            raise TenxRunError(f"No whitelist file available for {chem.id}")
        values: set[str] = set()
        with open_text_auto(path) as handle:
            for line in handle:
                fields = line.strip().split()
                if not fields:
                    continue
                bc = fields[0].upper()
                if len(bc) >= chem.cb_len:
                    values.add(bc[: chem.cb_len])
        if not values:
            raise TenxRunError(f"Whitelist {path} contains no >= {chem.cb_len} bp barcodes")
        self.cache[chem.id] = values
        return values


def match_whitelists(infos: List[FastqInfo], seq_samples: Dict[str, List[str]], store: WhitelistStore, args: argparse.Namespace) -> None:
    available = store.available()
    if not available:
        raise TenxRunError(f"No known 10x whitelist files found in {store.whitelist_dir}")
    for info in infos:
        seqs = seq_samples[info.path]
        for chem in available:
            wl = store.get(chem)
            for offset in chem.barcode_offsets:
                if info.common_length < offset + chem.cb_len:
                    continue
                hits = sum(1 for seq in seqs if len(seq) >= offset + chem.cb_len and seq[offset : offset + chem.cb_len] in wl)
                frac = hits / len(seqs)
                if hits >= args.min_whitelist_hits and frac >= args.min_whitelist_fraction:
                    info.matches.append(Match(chem.id, chem.modality, chem.label, offset, hits, len(seqs), frac))
        info.matches.sort(key=lambda m: (m.fraction, m.hits, -m.offset), reverse=True)


def matches_for(info: FastqInfo, modality: Optional[str] = None) -> List[Match]:
    return [m for m in info.matches if modality is None or m.modality == modality]


def best_match(info: FastqInfo, modality: Optional[str] = None) -> Optional[Match]:
    ms = matches_for(info, modality)
    return ms[0] if ms else None


def geometry_gex_candidate(info: FastqInfo, args: argparse.Namespace) -> bool:
    m = best_match(info, "gex")
    if not m:
        return False
    chem = CHEM_BY_ID[m.chemistry_id]
    if chem.split_v1:
        if info.common_length == chem.cb_len:
            return True
        return (info.common_length - chem.cb_len) in args.allow_v1_umi_lengths
    return chem.umi_len is not None and info.common_length >= chem.cb_len + chem.umi_len


def geometry_atac_candidate(info: FastqInfo) -> bool:
    return best_match(info, "atac") is not None and info.common_length in ATAC_BARCODE_LENGTHS


def require_constant(info: FastqInfo, role: str, min_fraction: float) -> None:
    if info.common_fraction < min_fraction:
        raise TenxRunError(f"{info.basename}: {role} has variable lengths {info.length_counts}; unsafe to rename")


# ---------------------------------------------------------------------------
# Role inference
# ---------------------------------------------------------------------------


def choose_one(items: List[FastqInfo], what: str) -> FastqInfo:
    if len(items) != 1:
        raise TenxRunError(f"Expected exactly one {what}, found {len(items)}: {[x.basename for x in items]}")
    return items[0]


def assign_index_roles(index_files: List[FastqInfo], warnings: List[str]) -> Dict[str, FastqInfo]:
    if not index_files:
        return {}
    if any(f.common_length not in INDEX_LENGTHS for f in index_files):
        raise TenxRunError(f"Non-index-like files remain: {[f'{x.basename}:{x.common_length}' for x in index_files]}")
    if len(index_files) > 2:
        raise TenxRunError(f"More than two index-like files remain: {[x.basename for x in index_files]}")
    roles: Dict[str, FastqInfo] = {}
    for f in index_files:
        if f.explicit_role in {"I1", "I2"}:
            if f.explicit_role in roles:
                raise TenxRunError(f"Duplicate {f.explicit_role} index hints among run FASTQs")
            roles[f.explicit_role] = f
    leftovers = [f for f in sorted(index_files, key=lambda x: natural_key(x.path)) if f not in roles.values()]
    if leftovers:
        if "I1" not in roles:
            roles["I1"] = leftovers.pop(0)
        if leftovers and "I2" not in roles:
            roles["I2"] = leftovers.pop(0)
        if len(index_files) == 2 and any(f.explicit_role not in {"I1", "I2"} for f in index_files):
            warnings.append("I1/I2 assignment used file order because one or both index files lacked explicit I1/I2 names")
    return roles


def classify_gex(infos: List[FastqInfo], store: WhitelistStore, args: argparse.Namespace, allow_excluded: bool) -> Tuple[ChemistryDef, Dict[str, FastqInfo], List[OutputEntry], List[str], List[str], Optional[int]]:
    warnings: List[str] = []
    excluded: List[str] = []
    barcode_candidates = [i for i in infos if geometry_gex_candidate(i, args)]
    barcode = choose_one(barcode_candidates, "GEX barcode read with valid CB+UMI geometry")
    match = best_match(barcode, "gex")
    assert match is not None
    chem = CHEM_BY_ID[match.chemistry_id]
    require_constant(barcode, "GEX barcode read", args.min_barcode_common_fraction)
    roles: Dict[str, FastqInfo] = {}
    outputs: List[OutputEntry] = []

    if chem.split_v1 and barcode.common_length == chem.cb_len:
        remaining = [i for i in infos if i is not barcode]
        umi_files = [i for i in remaining if i.common_length in args.allow_v1_umi_lengths]
        umi = choose_one(umi_files, "3' v1 UMI read")
        require_constant(umi, "3' v1 UMI read", args.min_barcode_common_fraction)
        bio_files = [i for i in remaining if i is not umi and i.common_length >= args.min_bio_read_length]
        bio = choose_one(bio_files, "3' v1 biological read")
        index_files = [i for i in remaining if i is not umi and i is not bio]
        index_roles = assign_index_roles(index_files, warnings)
        roles = {"CB": barcode, "UMI": umi, "R2": bio, **index_roles}
        roles["CB"].inferred_role = "CB"
        roles["UMI"].inferred_role = "UMI"
        roles["R2"].inferred_role = "R2"
        for rt, fi in index_roles.items():
            fi.inferred_role = rt
        umi_len = umi.common_length
        outputs.append(OutputEntry("R1", cr_name(args.fastq_prefix, "R1"), [barcode.path, umi.path], [barcode.common_length, umi.common_length], True, "synthetic CB+UMI for 3' v1"))
        outputs.append(OutputEntry("R2", cr_name(args.fastq_prefix, "R2"), [bio.path], [bio.common_length]))
        for rt in ("I1", "I2"):
            if rt in index_roles:
                outputs.append(OutputEntry(rt, cr_name(args.fastq_prefix, rt), [index_roles[rt].path], [index_roles[rt].common_length]))
        return chem, roles, outputs, warnings, excluded, umi_len

    if chem.split_v1:
        umi_len = barcode.common_length - chem.cb_len
        if umi_len not in args.allow_v1_umi_lengths:
            raise TenxRunError(f"Merged 3' v1 read length {barcode.common_length} implies UMI {umi_len}, not in {sorted(args.allow_v1_umi_lengths)}")
    else:
        expected = chem.cb_len + (chem.umi_len or 0)
        if barcode.common_length < expected:
            raise TenxRunError(f"{barcode.basename}: GEX barcode read length {barcode.common_length} is shorter than CB+UMI {expected}")
        umi_len = chem.umi_len
        if barcode.common_length > expected:
            warnings.append(f"Barcode read {barcode.basename} is over-sequenced: {barcode.common_length} bp > CB+UMI {expected} bp")

    remaining = [i for i in infos if i is not barcode]
    long_files = [i for i in remaining if i.common_length >= args.min_bio_read_length]
    if allow_excluded:
        preferred = [i for i in long_files if i.common_length >= args.min_gex_r2_length]
        if len(preferred) == 1:
            excluded.extend(i.path for i in long_files if i not in preferred)
            long_files = preferred
    bio = choose_one(long_files, "GEX biological R2 read")
    require_constant(bio, "GEX biological read", args.min_bio_common_fraction)
    index_files = [i for i in remaining if i is not bio and i.path not in excluded]
    if allow_excluded and any(geometry_atac_candidate(i) for i in infos):
        keep = [i for i in index_files if i.common_length == 10]
        excluded.extend(i.path for i in index_files if i not in keep)
        index_files = keep
    index_roles = assign_index_roles(index_files, warnings)
    roles = {"R1": barcode, "R2": bio, **index_roles}
    for rt, fi in roles.items():
        fi.inferred_role = rt
    outputs.append(OutputEntry("R1", cr_name(args.fastq_prefix, "R1"), [barcode.path], [barcode.common_length]))
    outputs.append(OutputEntry("R2", cr_name(args.fastq_prefix, "R2"), [bio.path], [bio.common_length]))
    for rt in ("I1", "I2"):
        if rt in index_roles:
            outputs.append(OutputEntry(rt, cr_name(args.fastq_prefix, rt), [index_roles[rt].path], [index_roles[rt].common_length]))
    return chem, roles, outputs, warnings, excluded, umi_len


def classify_atac(infos: List[FastqInfo], store: WhitelistStore, args: argparse.Namespace) -> Tuple[ChemistryDef, Dict[str, FastqInfo], List[OutputEntry], List[str], List[str], Optional[int]]:
    warnings: List[str] = []
    barcode_candidates = [i for i in infos if geometry_atac_candidate(i)]
    barcode = choose_one(barcode_candidates, "ATAC barcode-index read")
    match = best_match(barcode, "atac")
    assert match is not None
    chem = CHEM_BY_ID[match.chemistry_id]
    require_constant(barcode, "ATAC barcode-index read", args.min_barcode_common_fraction)
    if barcode.common_length == 24 and match.offset not in (0, 8):
        raise TenxRunError(f"{barcode.basename}: unexpected ATAC barcode offset {match.offset}")

    remaining = [i for i in infos if i is not barcode]
    genomic = [i for i in remaining if i.common_length >= args.min_atac_genomic_read_length]
    if len(genomic) != 2:
        raise TenxRunError(f"Expected two ATAC genomic reads, found {len(genomic)}: {[x.basename for x in genomic]}")
    explicit_r1 = [i for i in genomic if i.explicit_role == "R1"]
    explicit_r2 = [i for i in genomic if i.explicit_role in {"R2", "R3", "R4"}]
    if len(explicit_r1) == 1 and len(explicit_r2) == 1:
        r1, r2 = explicit_r1[0], explicit_r2[0]
    elif genomic[0].common_length != genomic[1].common_length:
        ordered = sorted(genomic, key=lambda x: x.common_length, reverse=True)
        r1, r2 = ordered[0], ordered[1]
    elif all(i.ordinal is not None for i in genomic):
        ordered = sorted(genomic, key=lambda x: x.ordinal or 0)
        r1, r2 = ordered[0], ordered[1]
        warnings.append("ATAC genomic reads had equal length; used numeric read order to assign R1/R2")
    else:
        raise TenxRunError("ATAC genomic reads have equal length and no usable R1/R2/R3/numeric ordering hint")
    leftover = [i for i in remaining if i is not r1 and i is not r2]
    index_roles = assign_index_roles(leftover, warnings)
    # ATAC cell barcode is emitted as I2; any ordinary sample index is I1.
    if "I2" in index_roles:
        raise TenxRunError("Found a separate I2 hint in addition to the ATAC barcode-index read")
    roles = {"R1": r1, "R2": r2, "I2": barcode, **index_roles}
    for rt, fi in roles.items():
        fi.inferred_role = rt
    outputs = [
        OutputEntry("R1", cr_name(args.fastq_prefix, "R1"), [r1.path], [r1.common_length]),
        OutputEntry("R2", cr_name(args.fastq_prefix, "R2"), [r2.path], [r2.common_length]),
        OutputEntry("I2", cr_name(args.fastq_prefix, "I2"), [barcode.path], [barcode.common_length], False, "ATAC cell barcode index"),
    ]
    if "I1" in index_roles:
        outputs.append(OutputEntry("I1", cr_name(args.fastq_prefix, "I1"), [index_roles["I1"].path], [index_roles["I1"].common_length]))
    return chem, roles, outputs, warnings, [], None


def cr_name(prefix: str, read_type: str) -> str:
    return f"{prefix}_{read_type}_001.fastq.gz"


def infer_plan(infos: List[FastqInfo], store: WhitelistStore, args: argparse.Namespace) -> Plan:
    has_gex = any(geometry_gex_candidate(i, args) for i in infos)
    has_atac = any(geometry_atac_candidate(i) for i in infos)
    if has_gex and has_atac and not args.prefer_gex:
        raise TenxRunError("Both GEX and ATAC barcode evidence detected in one run; refusing to mix modalities without --prefer-gex")
    if has_gex:
        chem, roles, outputs, warnings, excluded, umi_len = classify_gex(infos, store, args, allow_excluded=args.prefer_gex)
    elif has_atac:
        if args.accept_modalities == "gex":
            raise TenxRunError("Detected ATAC-only 10x run; rerun with --accept-modalities both/atac to audit it, but do not feed it to STARsolo GEX")
        chem, roles, outputs, warnings, excluded, umi_len = classify_atac(infos, store, args)
    else:
        best = []
        for i in infos:
            if i.matches:
                m = i.matches[0]
                best.append(f"{i.basename}:{m.chemistry_id}:{m.fraction:.1%}:len={i.common_length}")
        detail = "; ".join(best) or "no whitelist hits"
        raise TenxRunError(f"No supported 10x run layout matched whitelist evidence plus read geometry ({detail})")

    if args.accept_modalities != "both" and chem.modality != args.accept_modalities:
        raise TenxRunError(f"Detected modality {chem.modality}, but --accept-modalities={args.accept_modalities}")

    read_lengths: Dict[str, int] = {}
    for out in outputs:
        if out.synthetic:
            read_lengths[out.read_type] = sum(out.source_lengths)
        else:
            read_lengths[out.read_type] = out.source_lengths[0]
    role_map = {Path(fi.path).name: role for role, fi in roles.items()}
    confidence = "layout_only" if chem.id in {"gex_3pv2_or_5pv1v2", "gex_3pv3_family"} else "unique"
    return Plan(
        ok=True,
        run_id=args.run_id,
        sample_id=args.sample_id,
        fastq_prefix=args.fastq_prefix,
        modality=chem.modality,
        chemistry_id=chem.id,
        chemistry_label=chem.label,
        chemistry_version=chem.version,
        chemistry_group=chem.chemistry_group,
        whitelist=store.path_for(chem.id),
        cb_len=chem.cb_len,
        umi_len=umi_len,
        strand_hint=chem.strand_hint,
        layout_confidence=confidence,
        read_lengths=read_lengths,
        source_role_map=role_map,
        outputs=outputs,
        excluded_files=excluded,
        warnings=warnings,
        stats=infos,
        checks={},
    )


# ---------------------------------------------------------------------------
# Validation and output writing
# ---------------------------------------------------------------------------


def validate_plan(plan: Plan, args: argparse.Namespace) -> Dict[str, object]:
    checks: Dict[str, object] = {}
    selected_sources = []
    for out in plan.outputs:
        for src in out.sources:
            if src not in selected_sources:
                selected_sources.append(src)
    if args.check_counts:
        counts = {Path(src).name: count_records_strict(Path(src)) for src in selected_sources}
        if len(set(counts.values())) != 1:
            raise TenxRunError(f"Read counts differ across selected run FASTQs: {counts}")
        checks["record_counts"] = counts
    else:
        checks["record_counts"] = "skipped"
    if args.check_ids:
        ref = Path(selected_sources[0])
        ref_ids = first_ids(ref, args.id_check_records)
        mismatches: Dict[str, int] = {}
        for src in selected_sources[1:]:
            ids = first_ids(Path(src), args.id_check_records)
            n = min(len(ref_ids), len(ids))
            mm = sum(1 for a, b in zip(ref_ids[:n], ids[:n]) if a != b)
            if mm:
                mismatches[Path(src).name] = mm
        if mismatches:
            raise TenxRunError(f"Read IDs differ across selected run FASTQs: {mismatches}")
        checks["id_concordance_records"] = len(ref_ids)
    else:
        checks["id_concordance_records"] = "skipped"
    return checks


def place_as_gz(src: Path, dest: Path, action: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    comp = compression_of(src)
    if comp == "gz":
        if action == "symlink":
            os.symlink(os.path.abspath(src), dest)
            return
        if action == "hardlink":
            try:
                os.link(src.resolve(), dest)
                return
            except OSError:
                shutil.copy2(src.resolve(), dest)
                return
        shutil.copy2(src.resolve(), dest)
        return
    with open_bytes_auto(src) as inp, gzip.open(dest, "wb", compresslevel=6) as out:
        shutil.copyfileobj(inp, out, length=1024 * 1024)


def merge_cb_umi(cb: Path, umi: Path, dest: Path, cb_len: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cb_iter = iter_fastq_records(cb)
    umi_iter = iter_fastq_records(umi)
    n = 0
    with gzip.open(dest, "wt", encoding="utf-8") as out:
        while True:
            try:
                cb_rec = next(cb_iter)
            except StopIteration:
                try:
                    next(umi_iter)
                    raise TenxRunError("UMI FASTQ has more records than CB FASTQ")
                except StopIteration:
                    break
            try:
                umi_rec = next(umi_iter)
            except StopIteration:
                raise TenxRunError("CB FASTQ has more records than UMI FASTQ")
            h1, s1, p1, q1 = cb_rec
            h2, s2, _p2, q2 = umi_rec
            if normalize_read_id(h1) != normalize_read_id(h2):
                raise TenxRunError(f"CB/UMI record ID mismatch at record {n+1}: {h1} vs {h2}")
            out.write(f"{h1}\n{s1[:cb_len]}{s2}\n{p1}\n{q1[:cb_len]}{q2}\n")
            n += 1


def emit_outputs(plan: Plan, outdir: Path, action: str, dry_run: bool) -> None:
    if dry_run:
        return
    names = [o.output_name for o in plan.outputs]
    if len(names) != len(set(names)):
        raise TenxRunError(f"Duplicate output FASTQ names in run plan: {names}")
    for out in plan.outputs:
        dest = outdir / out.output_name
        if out.synthetic:
            if len(out.sources) != 2:
                raise TenxRunError(f"Synthetic output requires two source files: {out}")
            merge_cb_umi(Path(out.sources[0]), Path(out.sources[1]), dest, cb_len=plan.cb_len or 0)
        else:
            if len(out.sources) != 1:
                raise TenxRunError(f"Non-synthetic output requires one source file: {out}")
            place_as_gz(Path(out.sources[0]), dest, action)


def plan_to_json(plan: Optional[Plan], error: Optional[str] = None, stats: Optional[List[FastqInfo]] = None) -> Dict[str, object]:
    if plan is None:
        return {"ok": False, "error": error, "stats": [s.to_json() for s in (stats or [])], "tool": {"name": "infer_10x_run.py"}}
    return {
        "ok": True,
        "run_id": plan.run_id,
        "sample_id": plan.sample_id,
        "fastq_prefix": plan.fastq_prefix,
        "modality": plan.modality,
        "chemistry_id": plan.chemistry_id,
        "chemistry_label": plan.chemistry_label,
        "chemistry_version": plan.chemistry_version,
        "chemistry_group": plan.chemistry_group,
        "whitelist": plan.whitelist,
        "cb_len": plan.cb_len,
        "umi_len": plan.umi_len,
        "strand_hint": plan.strand_hint,
        "layout_confidence": plan.layout_confidence,
        "read_lengths": plan.read_lengths,
        "source_role_map": plan.source_role_map,
        "outputs": [asdict(o) for o in plan.outputs],
        "excluded_files": plan.excluded_files,
        "warnings": plan.warnings,
        "checks": plan.checks,
        "stats": [s.to_json() for s in plan.stats],
        "tool": {"name": "infer_10x_run.py"},
    }


def write_manifest_tsv(plan: Plan, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        out.write("run_id\tsample_id\tfastq_prefix\tmodality\tchemistry_id\tread_type\tread_length\toutput_name\tsynthetic\tsources\n")
        for e in plan.outputs:
            out.write(
                f"{plan.run_id}\t{plan.sample_id}\t{plan.fastq_prefix}\t{plan.modality}\t{plan.chemistry_id}\t{e.read_type}\t"
                f"{plan.read_lengths[e.read_type]}\t{e.output_name}\t{e.synthetic}\t{';'.join(e.sources)}\n"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def derive_prefix(args: argparse.Namespace) -> str:
    if args.fastq_prefix:
        prefix = args.fastq_prefix
    else:
        sample = args.sample_id or args.run_id
        prefix = f"{sample}_S{args.sample_index}_L{args.lane:03d}"
    if not CR_PREFIX_RE.match(prefix):
        raise TenxRunError(f"FASTQ prefix must look like SAMPLE_S1_L001; got {prefix!r}")
    if not args.allow_unsafe_prefix and not SAFE_PREFIX_RE.match(prefix):
        raise TenxRunError(f"FASTQ prefix contains unsafe characters for Cell Ranger-style names: {prefix!r}")
    return prefix


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Infer one 10x run layout and emit Cell Ranger-style FASTQs")
    p.add_argument("--fastqs", nargs="+", type=Path, required=True, help="2-4 FASTQ/FASTQ.GZ/FASTQ.BZ2 files for one run")
    p.add_argument("--run-id", required=True, help="Source run id, e.g. SRR/ERR/run-lane id")
    p.add_argument("--sample-id", default=None, help="Sample id used if --fastq-prefix is not provided")
    p.add_argument("--fastq-prefix", default=None, help="CR prefix SAMPLE_S1_L001; usually computed upstream from run/lane/tag metadata")
    p.add_argument("--sample-index", type=int, default=1, help="S index used when deriving --fastq-prefix")
    p.add_argument("--lane", type=int, default=1, help="Lane number used when deriving --fastq-prefix")
    p.add_argument("--whitelist-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--json", dest="json_path", type=Path, default=Path("run_chemistry.json"))
    p.add_argument("--tsv", dest="tsv_path", type=Path, default=Path("run_chemistry.tsv"))
    p.add_argument("--accept-modalities", choices=("gex", "atac", "both"), default="gex")
    p.add_argument("--prefer-gex", action="store_true", help="For a wrongly mixed GEX+ATAC invocation, emit only the separable GEX run")
    p.add_argument("--sample-records", type=int, default=200000)
    p.add_argument("--min-whitelist-hits", type=int, default=1000)
    p.add_argument("--min-whitelist-fraction", type=float, default=0.20)
    p.add_argument("--min-bio-read-length", type=int, default=40)
    p.add_argument("--min-gex-r2-length", type=int, default=60)
    p.add_argument("--min-atac-genomic-read-length", type=int, default=40)
    p.add_argument("--min-barcode-common-fraction", type=float, default=0.98)
    p.add_argument("--min-bio-common-fraction", type=float, default=0.90)
    p.add_argument("--allow-v1-umi-lengths", default="5,10")
    p.add_argument("--action", choices=("copy", "symlink", "hardlink"), default="hardlink")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--no-check-counts", dest="check_counts", action="store_false", default=True)
    p.add_argument("--no-check-ids", dest="check_ids", action="store_false", default=True)
    p.add_argument("--id-check-records", type=int, default=1000)
    p.add_argument("--allow-unsafe-prefix", action="store_true")
    args = p.parse_args(argv)
    if not (2 <= len(args.fastqs) <= 4):
        p.error("This run-level script expects 2-4 FASTQ files")
    for f in args.fastqs:
        if not f.exists():
            p.error(f"FASTQ not found: {f}")
    if not args.whitelist_dir.exists():
        p.error(f"Whitelist directory not found: {args.whitelist_dir}")
    args.allow_v1_umi_lengths = {int(x) for x in str(args.allow_v1_umi_lengths).split(",") if x.strip()}
    if not args.sample_id:
        args.sample_id = args.run_id
    args.fastq_prefix = derive_prefix(args)
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    stats: List[FastqInfo] = []
    try:
        log("INFO", f"Run {args.run_id}: sampling {len(args.fastqs)} FASTQ files")
        stats, samples = collect_fastq_info(args.fastqs, args.sample_records)
        for s in stats:
            log("INFO", f"{s.basename}: len={s.common_length} ({s.common_fraction:.1%}), explicit={s.explicit_role}, ordinal={s.ordinal}, comp={s.compression}")
        store = WhitelistStore(args.whitelist_dir)
        if store.missing():
            log("WARN", "Missing optional whitelist definitions: " + ",".join(store.missing()))
        match_whitelists(stats, samples, store, args)
        plan = infer_plan(stats, store, args)
        plan.checks = validate_plan(plan, args)
        emit_outputs(plan, args.outdir, args.action, args.dry_run)
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json_path, "w", encoding="utf-8") as fh:
            json.dump(plan_to_json(plan), fh, indent=2, sort_keys=True)
        write_manifest_tsv(plan, args.tsv_path)
        for w in plan.warnings:
            log("WARN", w)
        if plan.excluded_files:
            log("WARN", "Excluded files: " + ", ".join(Path(x).name for x in plan.excluded_files))
        log("INFO", f"Detected {plan.chemistry_id} ({plan.modality}); wrote {len(plan.outputs)} output role(s)")
        return 0
    except Exception as exc:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json_path, "w", encoding="utf-8") as fh:
            json.dump(plan_to_json(None, error=str(exc), stats=stats), fh, indent=2, sort_keys=True)
        log("ERROR", str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
