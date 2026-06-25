#!/usr/bin/env python3
"""Infer 10x FASTQ read roles and rename files to Cell Ranger FASTQ names.

The detector is deliberately conservative: it combines observed read lengths,
read-id compatibility, and exact whitelist matches for the cell-barcode segment.
It refuses to rename when multiple candidate chemistries imply different read-role
mappings. Chemistry names that share the same FASTQ geometry are reported as
ambiguous while still allowing a safe Cell Ranger/STARsolo rename.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import re
import shutil
import statistics
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

READ_TYPES = ("R1", "R2", "I1", "I2")
MIN_RNA_BASES = 25

ALIASES: Mapping[str, Tuple[str, ...]] = {
    "737K-april-2014_rc": ("737K-april-2014_rc.txt", "737K-april-2014_rc"),
    "737K-august-2016": ("737K-august-2016.txt", "737K-august-2016"),
    "3M-february-2018": ("3M-february-2018.txt", "3M-february-2018"),
    "3M-february-2018_TRU": ("3M-february-2018.txt", "3M-february-2018", "3M-february-2018_TRU"),
    "3M-february-2018_NXT": ("3M-february-2018.txt", "3M-february-2018", "3M-february-2018_NXT"),
    "3M-3pgex-may-2023": ("3M-3pgex-may-2023.txt", "3M-3pgex-may-2023"),
    "3M-3pgex-may-2023_TRU": ("3M-3pgex-may-2023.txt", "3M-3pgex-may-2023", "3M-3pgex-may-2023_TRU"),
    "3M-3pgex-may-2023_NXT": ("3M-3pgex-may-2023.txt", "3M-3pgex-may-2023", "3M-3pgex-may-2023_NXT"),
    "3M-5pgex-jan-2023": ("3M-5pgex-jan-2023.txt", "3M-5pgex-jan-2023"),
    "737K-arc-v1": ("737K-arc-v1.txt", "737K-arc-v1"),
    "9K-LT-march-2021": ("9K-LT-march-2021.txt", "9K-LT-march-2021"),
    "737K-fixed-rna-profiling": ("737K-fixed-rna-profiling.txt", "737K-fixed-rna-profiling"),
    "737K-flex-v2": ("737K-flex-v2.txt", "737K-flex-v2"),
}


def chem(name: str, desc: str, wl: str, bc_rt: str, bc_len: int, umi_rt: str, umi_off: int, umi_len: int,
         rna_rt: str, rna_off: int = 0, rna_min: int = MIN_RNA_BASES,
         rna2_rt: Optional[str] = None, rna2_off: int = 0, rna2_min: int = MIN_RNA_BASES,
         umi_min: Optional[int] = None, star: bool = True) -> Dict[str, object]:
    return dict(name=name, description=desc, whitelist=wl, bc_read=bc_rt, bc_len=bc_len, bc_off=0,
                umi_read=umi_rt, umi_off=umi_off, umi_len=umi_len, umi_min=umi_len if umi_min is None else umi_min,
                rna_read=rna_rt, rna_off=rna_off, rna_min=rna_min,
                rna2_read=rna2_rt, rna2_off=rna2_off, rna2_min=rna2_min, starsolo=star)


def catalog() -> List[Dict[str, object]]:
    c: List[Dict[str, object]] = []
    # 3' GEX and Multiome GEX
    c += [
        chem("SC3Pv1", "Single Cell 3' v1", "737K-april-2014_rc", "I1", 14, "R2", 0, 10, "R1"),
        chem("SC3Pv1-reconstructed", "Single Cell 3' v1 reconstructed CB+UMI", "737K-april-2014_rc", "R1", 14, "R1", 14, 10, "R2"),
        chem("SC3Pv2", "Single Cell 3' v2", "737K-august-2016", "R1", 16, "R1", 16, 10, "R2"),
    ]
    for suffix, wl in (("polyA", "3M-february-2018_TRU"), ("CS1", "3M-february-2018_NXT")):
        for extra in ("", "-OCM", "HT"):
            n = f"SC3Pv3{extra}-{suffix}" if extra == "HT" else f"SC3Pv3-{suffix}{extra}"
            d = f"Single Cell 3' v3 {extra or ''} ({suffix})".replace("  ", " ")
            c.append(chem(n, d, wl, "R1", 16, "R1", 16, 12, "R2", umi_min=10))
    c.append(chem("SC3Pv3LT", "Single Cell 3' v3 LT", "9K-LT-march-2021", "R1", 16, "R1", 16, 12, "R2", umi_min=10))
    for suffix, wl in (("polyA", "3M-3pgex-may-2023_TRU"), ("CS1", "3M-3pgex-may-2023_NXT")):
        c.append(chem(f"SC3Pv4-{suffix}", f"Single Cell 3' v4 ({suffix})", wl, "R1", 16, "R1", 16, 12, "R2", umi_min=10))
        c.append(chem(f"SC3Pv4-{suffix}-OCM", f"Single Cell 3' v4 OCM ({suffix})", wl, "R1", 16, "R1", 16, 12, "R2", umi_min=10))
    c.append(chem("ARC-v1", "Single Cell Multiome ATAC + Gene Expression v1", "737K-arc-v1", "R1", 16, "R1", 16, 12, "R2", umi_min=10))

    # 5' GEX and V(D)J. R1-only and VDJ are included so unsafe STARsolo/simple GEX cases are visible in the report.
    for n, d, wl, ulen, r1off in [
        ("SC5P-R1", "Single Cell 5' R1-only", "737K-august-2016", 10, 41),
        ("SC5P-R1-v3", "Single Cell 5' R1-only v3", "3M-5pgex-jan-2023", 12, 43),
        ("SC5P-R1-OCM-v3", "Single Cell 5' R1-only OCM v3", "3M-5pgex-jan-2023", 12, 43),
    ]:
        c.append(chem(n, d, wl, "R1", 16, "R1", 16, ulen, "R1", rna_off=r1off, umi_min=10, star=False))
    for n, d, wl, ulen in [
        ("SC5P-R2", "Single Cell 5' R2-only", "737K-august-2016", 10),
        ("SC5P-R2-OCM", "Single Cell 5' R2-only OCM", "737K-august-2016", 10),
        ("SC5PHT", "Single Cell 5' HT", "737K-august-2016", 10),
        ("SC5P-R2-v3", "Single Cell 5' R2-only v3", "3M-5pgex-jan-2023", 12),
        ("SC5P-R2-OCM-v3", "Single Cell 5' R2-only OCM v3", "3M-5pgex-jan-2023", 12),
        ("SCVDJ-R2", "Single Cell V(D)J R2-only", "737K-august-2016", 10),
        ("SCVDJ-R2-v3", "Single Cell V(D)J R2-only v3", "3M-5pgex-jan-2023", 12),
        ("SCVDJ-R2-OCM-v3", "Single Cell V(D)J R2-only OCM v3", "3M-5pgex-jan-2023", 12),
        ("SC-FB", "Single Cell 3' v2 or 5' Feature Barcode", "737K-august-2016", 10),
    ]:
        c.append(chem(n, d, wl, "R1", 16, "R1", 16, ulen, "R2", umi_min=10, star=not n.startswith("SCVDJ")))
    for n, d, wl, ulen, r1min in [
        ("SC5P-PE", "Single Cell 5' PE", "737K-august-2016", 10, 81),
        ("SC5P-PE-v3", "Single Cell 5' PE v3", "3M-5pgex-jan-2023", 12, 83),
        ("SC5P-PE-OCM-v3", "Single Cell 5' PE OCM v3", "3M-5pgex-jan-2023", 12, 83),
        ("SCVDJ", "Single Cell V(D)J", "737K-august-2016", 10, 66),
        ("SCVDJ-v3", "Single Cell V(D)J v3", "3M-5pgex-jan-2023", 12, 68),
        ("SCVDJ-v3-OCM", "Single Cell V(D)J OCM v3", "3M-5pgex-jan-2023", 12, 68),
    ]:
        c.append(chem(n, d, wl, "R1", 16, "R1", 16, ulen, "R1", rna_off=0, rna_min=r1min, rna2_rt="R2", umi_min=10, star=not n.startswith("SCVDJ")))

    # Fixed RNA Profiling / Flex. These use the same R1/R2 naming but are flagged as not simple STARsolo GEX.
    flex_specs = [
        ("SFRP", "Flex Gene Expression", "737K-fixed-rna-profiling", 30),
        ("SFRP-no-trim-R2", "SFRP with untrimmed R2", "737K-fixed-rna-profiling", MIN_RNA_BASES),
        ("MFRP-RNA", "Flex Gene Expression", "737K-fixed-rna-profiling", 50),
        ("MFRP-Ab", "Flex Antibody", "737K-fixed-rna-profiling", 50),
        ("MFRP-RNA-R1", "Flex Gene Expression probe barcode on R1", "737K-fixed-rna-profiling", 30),
        ("MFRP-Ab-R1", "Flex Antibody probe barcode on R1", "737K-fixed-rna-profiling", 30),
        ("MFRP-R1-48-uncollapsed", "Fixed RNA profiling probe barcode on R1", "737K-fixed-rna-profiling", 30),
        ("MFRP-47", "Fixed RNA profiling 47 probe barcodes", "737K-fixed-rna-profiling", 50),
        ("MFRP-96-R1", "Flex 96-plex beta", "737K-fixed-rna-profiling", 30),
        ("MFRP-96-RNA-R2", "Flex Gene Expression 96-plex beta", "737K-fixed-rna-profiling", MIN_RNA_BASES),
        ("MFRP-uncollapsed", "Multiplex fixed RNA profiling", "737K-fixed-rna-profiling", 50),
        ("MFRP-Ab-R2pos50", "Flex Antibody probe barcode at R2:50", "737K-fixed-rna-profiling", 50),
        ("MFRP-CRISPR", "Flex CRISPR", "737K-fixed-rna-profiling", MIN_RNA_BASES),
        ("MFRP-R1-no-trim-R2", "MFRP probe barcode on R1 with untrimmed R2", "737K-fixed-rna-profiling", MIN_RNA_BASES),
        ("Flex-v2-singleplex", "GEM-X Flex v2 singleplex", "737K-flex-v2", 30),
        ("Flex-v2-96-R1", "Flex v2 96-plex beta", "737K-flex-v2", 30),
        ("Flex-v2-96-RNA-R2", "Flex v2 Gene Expression 96-plex beta", "737K-flex-v2", MIN_RNA_BASES),
        ("Flex-v2-R1", "GEM-X Flex v2", "737K-flex-v2", 30),
        ("Flex-v2-RNA-R2", "GEM-X Flex v2 Gene Expression", "737K-flex-v2", MIN_RNA_BASES),
        ("Flex-v2-Ab-R2:45", "GEM-X Flex v2 Antibody R2:45", "737K-flex-v2", MIN_RNA_BASES),
        ("Flex-v2-Ab-R2:64", "GEM-X Flex v2 Antibody R2:64", "737K-flex-v2", MIN_RNA_BASES),
        ("Flex-v2-CRISPR-R2:1", "GEM-X Flex v2 CRISPR R2:1", "737K-flex-v2", MIN_RNA_BASES),
    ]
    for n, d, wl, r2min in flex_specs:
        c.append(chem(n, d, wl, "R1", 16, "R1", 16, 12, "R2", rna_min=r2min, umi_min=10, star=False))
    return c


class DetectionError(RuntimeError):
    pass


class Whitelists:
    def __init__(self, root: Path):
        self.root = root
        self.cache: Dict[str, Optional[set[str]]] = {}
        self.paths: Dict[str, str] = {}

    def load(self, name: str) -> Optional[set[str]]:
        if name in self.cache:
            return self.cache[name]
        for candidate in ALIASES.get(name, (f"{name}.txt", name)):
            path = self.root / candidate
            if path.exists():
                vals = {line.split()[0] for line in path.read_text().splitlines() if line.strip()}
                self.cache[name] = vals
                self.paths[name] = str(path)
                return vals
        self.cache[name] = None
        return None


def wl_key(offset: int, length: int, whitelist: str) -> str:
    return f"{offset}:{length}:{whitelist}"


def openfq(path: Path) -> Iterator[str]:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            yield from (x.rstrip("\n") for x in handle)
    else:
        with path.open() as handle:
            yield from (x.rstrip("\n") for x in handle)


def parse_name(path: Path) -> Tuple[Optional[str], str, Optional[int], Optional[str]]:
    name = path.name
    m = re.match(r"^(.+?)_S\d+(?:_L(\d{3}))?_([RI][12])_\d{3}\.f(?:ast)?q(?:\.gz)?$", name, re.I)
    if m:
        lane = int(m.group(2) or 1)
        return m.group(1), f"{m.group(1)}_L{lane:03d}", lane, m.group(3).upper()
    m = re.match(r"^(.+?)[._-](R[123]|I[12]|[1234])(?:[._-].*)?\.f(?:ast)?q(?:\.gz)?$", name, re.I)
    if m:
        raw = m.group(2).upper()
        read = {"1": "R1", "2": "R2", "3": "R3", "4": "R4"}.get(raw, raw)
        return None, m.group(1), None, read
    stem = re.sub(r"\.f(?:ast)?q(?:\.gz)?$", "", name, flags=re.I)
    return None, stem, None, None


def header_core(header: str) -> str:
    x = header[1:] if header.startswith("@") else header
    x = x.split()[0]
    return re.sub(r"([/._-][1234])$", "", x)


def stats_for(path: Path, wl: Whitelists, specs: Sequence[Tuple[int, int, str]], sample_reads: int) -> Dict[str, object]:
    sample, group, lane, read = parse_name(path)
    length_counts: Dict[int, int] = collections.Counter()
    hits: Dict[str, int] = collections.Counter()
    heads: List[str] = []
    loaded = []
    for off, ln, wname in specs:
        vals = wl.load(wname)
        if vals:
            loaded.append((off, ln, wname, vals, wl_key(off, ln, wname)))
    sampled = 0
    rec = 0
    current_header = ""
    for line in openfq(path):
        if rec == 0:
            current_header = line
        elif rec == 1:
            seq = line.strip()
            sampled += 1
            length_counts[len(seq)] += 1
            if len(heads) < 256:
                heads.append(current_header)
            for off, ln, wname, vals, key in loaded:
                if len(seq) >= off + ln and seq[off:off + ln] in vals:
                    hits[key] += 1
            if sampled >= sample_reads:
                break
        rec = (rec + 1) % 4
    if sampled == 0:
        raise DetectionError(f"{path}: no FASTQ records")
    lengths = [x for ln, n in length_counts.items() for x in [ln] * n]
    return dict(path=path, sample=sample, group=group, lane=lane, name_read=read, n=sampled,
                lengths=dict(length_counts), median=int(statistics.median(lengths)),
                mode=max(length_counts.items(), key=lambda kv: (kv[1], kv[0]))[0],
                variable=len(length_counts) > 1, hits=dict(hits), headers=heads)


def headers_match(a: Mapping[str, object], b: Mapping[str, object]) -> bool:
    ah, bh = a["headers"], b["headers"]
    n = min(len(ah), len(bh), 100)
    if n == 0:
        return True
    return sum(1 for i in range(n) if header_core(ah[i]) == header_core(bh[i])) / n >= 0.95


def priority(c: Mapping[str, object]) -> Tuple[int, str]:
    n = str(c["name"])
    penalty = 0
    if "OCM" in n: penalty += 10
    if "CS1" in n: penalty += 1
    if n.startswith("SC-FB"): penalty += 20
    if n.startswith("SCVDJ"): penalty += 30
    if n.startswith(("SFRP", "MFRP", "Flex")): penalty += 40
    if n.endswith("reconstructed"): penalty += 50
    return penalty, n


def min_len(c: Mapping[str, object], read: str) -> int:
    if c["name"] == "SC5P-PE" and read == "R1":
        return 81
    if c["name"] in {"SC5P-PE-v3", "SC5P-PE-OCM-v3"} and read == "R1":
        return 83
    m = 0
    if c["bc_read"] == read:
        m = max(m, int(c["bc_off"]) + int(c["bc_len"]))
    if c["umi_read"] == read:
        m = max(m, int(c["umi_off"]) + int(c["umi_min"]))
    if c["rna_read"] == read:
        m = max(m, int(c["rna_off"]) + int(c["rna_min"]))
    if c.get("rna2_read") == read:
        m = max(m, int(c.get("rna2_off", 0)) + int(c.get("rna2_min", MIN_RNA_BASES)))
    return m


def choose_rna(stats: Sequence[Mapping[str, object]], required: int, used: Sequence[Path]) -> Optional[Mapping[str, object]]:
    used_set = set(used)
    possible = [s for s in stats if s["path"] not in used_set and int(s["median"]) >= required]
    return sorted(possible, key=lambda s: (int(s["median"]), str(s["path"])), reverse=True)[0] if possible else None


def call_one(c: Mapping[str, object], stats: Sequence[Mapping[str, object]], min_frac: float, min_hits: int, target: str, strict_names: bool) -> Optional[Dict[str, object]]:
    key = wl_key(0, int(c["bc_len"]), str(c["whitelist"]))
    bc_files = []
    for s in stats:
        h = int(s["hits"].get(key, 0))
        if h >= min_hits and h / int(s["n"]) >= min_frac and int(s["median"]) >= int(c["bc_len"]):
            bc_files.append(s)
    calls = []
    for bc in sorted(bc_files, key=lambda s: int(s["hits"].get(key, 0)), reverse=True):
        if strict_names and bc["name_read"] and bc["name_read"] != c["bc_read"]:
            continue
        roles: Dict[str, Path] = {str(c["bc_read"]): bc["path"]}
        used: List[Path] = [bc["path"]]
        ok = True
        if c["umi_read"] == c["bc_read"]:
            if int(bc["median"]) < int(c["umi_off"]) + int(c["umi_min"]):
                ok = False
        else:
            opts = [s for s in stats if s["path"] not in used and not s["variable"] and int(s["mode"]) == int(c["umi_len"])]
            if strict_names:
                opts = [s for s in opts if not s["name_read"] or s["name_read"] == c["umi_read"] or (c["name"] == "SC3Pv1" and s["name_read"] in {"R2", "R3"})]
            if not opts:
                ok = False
            else:
                chosen = opts[0]
                roles[str(c["umi_read"])] = chosen["path"]
                used.append(chosen["path"])
        if not ok:
            continue
        for rt in (str(c["rna_read"]), c.get("rna2_read")):
            if not rt:
                continue
            required = min_len(c, rt)
            if rt in roles:
                if int(next(s for s in stats if s["path"] == roles[rt])["median"]) < required:
                    ok = False
                    break
            else:
                chosen = choose_rna(stats, required, used)
                if not chosen or (strict_names and chosen["name_read"] and chosen["name_read"] != rt):
                    ok = False
                    break
                roles[rt] = chosen["path"]
                used.append(chosen["path"])
        if not ok:
            continue
        if not all(headers_match(bc, next(s for s in stats if s["path"] == p)) for p in used[1:]):
            continue
        synthetic = target == "starsolo" and c["name"] == "SC3Pv1"
        mapped = dict(roles)
        if synthetic:
            if not {"I1", "R2", "R1"}.issubset(mapped):
                continue
            mapped = {"R1_CB": mapped["I1"], "R1_UMI": mapped["R2"], "R2": mapped["R1"]}
        calls.append(dict(chem=c, roles=mapped, used=tuple(used), synthetic=synthetic,
                          score=int(bc["hits"].get(key, 0)) / int(bc["n"])))
    return sorted(calls, key=lambda x: (-x["score"], priority(x["chem"])))[0] if calls else None


def infer_group(stats: Sequence[Mapping[str, object]], chems: Sequence[Mapping[str, object]], args) -> Dict[str, object]:
    calls = [x for c in chems if (x := call_one(c, stats, args.min_whitelist_fraction, args.min_whitelist_hits, args.target, args.strict_names))]
    if not calls:
        detail = "; ".join(f"{Path(s['path']).name}: median={s['median']} read={s['name_read']} hits={s['hits']}" for s in stats)
        raise DetectionError("No 10x chemistry matched FASTQs. " + detail)
    sigs = {tuple(sorted((k, str(v)) for k, v in call["roles"].items())) + (("synthetic", str(call["synthetic"])),) for call in calls}
    if len(sigs) != 1:
        raise DetectionError("Ambiguous read-role mapping; refusing to rename. Candidates: " + "; ".join(str(c["chem"]["name"]) for c in calls))
    calls.sort(key=lambda x: (-x["score"], priority(x["chem"])))
    names = sorted({str(x["chem"]["name"]) for x in calls})
    warnings = []
    if len(names) > 1:
        msg = "Chemistry is not unique from FASTQ geometry alone; read-role mapping is unique. Candidates: " + ", ".join(names)
        if args.require_unique_chemistry:
            raise DetectionError(msg)
        warnings.append(msg)
    if args.target == "starsolo" and any(not x["chem"].get("starsolo", True) for x in calls):
        warnings.append("One or more matching chemistries are not simple STARsolo GEX layouts; confirm library type before alignment.")
    return dict(selected=calls[0], candidates=calls, warnings=warnings, ambiguous=len(names) > 1)


def records(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    it = iter(openfq(path))
    while True:
        try:
            yield next(it), next(it), next(it), next(it)
        except StopIteration:
            return


def link_or_copy(src: Path, dest: Path, mode: str) -> None:
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if not str(src).endswith(".gz"):
        with open(src, "rb") as inp, gzip.open(dest, "wb") as out:
            shutil.copyfileobj(inp, out)
    elif mode == "copy":
        shutil.copy2(src, dest)
    elif mode == "hardlink":
        try: os.link(src, dest)
        except OSError: os.symlink(src.resolve(), dest)
    else:
        os.symlink(src.resolve(), dest)


def merge_v1(cb: Path, umi: Path, out: Path) -> None:
    with gzip.open(out, "wt") as h:
        for i, (a, b) in enumerate(zip(records(cb), records(umi)), 1):
            if header_core(a[0]) != header_core(b[0]):
                raise DetectionError(f"v1 CB/UMI read id mismatch at record {i}: {a[0]} vs {b[0]}")
            h.write(f"{a[0]}\n{a[1]}{b[1]}\n{a[2]}\n{a[3]}{b[3]}\n")


def split_groups(stats: Sequence[Mapping[str, object]]) -> List[Tuple[str, List[Mapping[str, object]]]]:
    if all(s["name_read"] for s in stats):
        grouped: Dict[str, List[Mapping[str, object]]] = collections.defaultdict(list)
        for s in stats:
            grouped[str(s["group"])].append(s)
        return sorted(grouped.items())
    return [("all", list(stats))]


def out_name(sample: str, lane: int, rt: str, outdir: Path) -> Path:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample).strip("._-") or "sample"
    return outdir / f"{clean}_S1_L{lane:03d}_{rt}_001.fastq.gz"


def write_outputs(sample: str, lane: int, outdir: Path, group_stats: Sequence[Mapping[str, object]], result: Mapping[str, object], mode: str, keep_index: bool) -> List[Dict[str, str]]:
    outdir.mkdir(parents=True, exist_ok=True)
    call = result["selected"]
    roles: Mapping[str, Path] = call["roles"]
    renamed: List[Dict[str, str]] = []
    if call["synthetic"]:
        r1, r2 = out_name(sample, lane, "R1", outdir), out_name(sample, lane, "R2", outdir)
        merge_v1(Path(roles["R1_CB"]), Path(roles["R1_UMI"]), r1)
        link_or_copy(Path(roles["R2"]), r2, mode)
        renamed += [{"read_type":"R1", "source":f"{roles['R1_CB']}+{roles['R1_UMI']}", "dest":str(r1), "synthetic":"true"},
                    {"read_type":"R2", "source":str(roles["R2"]), "dest":str(r2), "synthetic":"false"}]
    else:
        for rt in READ_TYPES:
            if rt in roles:
                dest = out_name(sample, lane, rt, outdir)
                link_or_copy(Path(roles[rt]), dest, mode)
                renamed.append({"read_type":rt, "source":str(roles[rt]), "dest":str(dest), "synthetic":"false"})
        if keep_index:
            used = set(call["used"])
            done = {x["read_type"] for x in renamed}
            for s in group_stats:
                if s["path"] in used or s["name_read"] not in {"I1", "I2"} or s["name_read"] in done:
                    continue
                dest = out_name(sample, lane, s["name_read"], outdir)
                link_or_copy(Path(s["path"]), dest, mode)
                renamed.append({"read_type":s["name_read"], "source":str(s["path"]), "dest":str(dest), "synthetic":"false"})
    return renamed


def stat_json(s: Mapping[str, object]) -> Dict[str, object]:
    return {"path": str(s["path"]), "name_read_type": s["name_read"], "name_group": s["group"],
            "median_length": s["median"], "mode_length": s["mode"], "length_counts": s["lengths"],
            "sampled_reads": s["n"], "variable_lengths": s["variable"], "whitelist_hits": s["hits"]}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fastqs", nargs="+", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--whitelist-dir", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--target", choices=("starsolo", "cellranger"), default="starsolo")
    p.add_argument("--sample-reads", type=int, default=200000)
    p.add_argument("--min-whitelist-fraction", type=float, default=0.20)
    p.add_argument("--min-whitelist-hits", type=int, default=100)
    p.add_argument("--copy-mode", choices=("symlink", "hardlink", "copy"), default="symlink")
    p.add_argument("--strict-names", action="store_true")
    p.add_argument("--require-unique-chemistry", action="store_true")
    p.add_argument("--no-index-reads", action="store_true")
    p.add_argument("--json-report", type=Path, default=Path("10x_fastq_inference.json"))
    p.add_argument("--tsv-report", type=Path, default=Path("10x_fastq_inference.tsv"))
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    chems = catalog()
    specs = sorted({(0, int(c["bc_len"]), str(c["whitelist"])) for c in chems})
    wl = Whitelists(args.whitelist_dir)
    missing = [str(x) for x in args.fastqs if not x.exists()]
    if missing:
        raise DetectionError("Missing FASTQ inputs: " + ", ".join(missing))
    stats = [stats_for(x, wl, specs, args.sample_reads) for x in args.fastqs]
    groups = split_groups(stats)
    group_reports = []
    all_renamed: List[Dict[str, str]] = []
    for lane, (gname, gstats) in enumerate(groups, 1):
        res = infer_group(gstats, chems, args)
        renamed = write_outputs(args.sample, lane, args.outdir, gstats, res, args.copy_mode, not args.no_index_reads)
        all_renamed += renamed
        group_reports.append({
            "group": gname, "output_lane": lane,
            "selected_chemistry": res["selected"]["chem"]["name"],
            "selected_description": res["selected"]["chem"]["description"],
            "candidate_chemistries": sorted({x["chem"]["name"] for x in res["candidates"]}),
            "ambiguous_chemistry": res["ambiguous"], "safe_rename": True,
            "warnings": res["warnings"], "input_files": [stat_json(s) for s in gstats],
            "renamed_files": renamed,
        })
    common = set(group_reports[0]["candidate_chemistries"])
    for r in group_reports[1:]:
        common &= set(r["candidate_chemistries"])
    if not common:
        raise DetectionError("Multiple groups did not share any compatible 10x chemistry")
    by_name = {c["name"]: c for c in chems}
    selected = sorted(common, key=lambda n: priority(by_name[n]))[0]
    report = group_reports[0] if len(group_reports) == 1 else {
        "groups": group_reports,
        "input_files": [x for r in group_reports for x in r["input_files"]],
        "renamed_files": all_renamed,
        "candidate_chemistries": sorted({x for r in group_reports for x in r["candidate_chemistries"]}),
        "ambiguous_chemistry": any(r["ambiguous_chemistry"] for r in group_reports),
        "safe_rename": True,
        "warnings": sorted({w for r in group_reports for w in r["warnings"]}),
    }
    report.update({"sample": args.sample, "target": args.target, "selected_chemistry": selected,
                   "resolved_whitelists": wl.paths})
    args.json_report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.tsv_report.write_text("sample\tselected_chemistry\tambiguous_chemistry\tsafe_rename\twarnings\n" +
                               f"{args.sample}\t{selected}\t{report['ambiguous_chemistry']}\ttrue\t" + "; ".join(report.get("warnings", [])) + "\n")
    for w in report.get("warnings", []):
        print("WARNING:", w, file=sys.stderr)
    print(f"Selected chemistry: {selected}", file=sys.stderr)
    for x in all_renamed:
        print(f"{x['source']} -> {x['dest']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except DetectionError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(2)
