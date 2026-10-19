#!/usr/bin/env python3
"""
Collate one run's per-dataset, per-run and per-sample artefacts into run-level tables.

Replaces the collectFile blocks main.nf used to carry. Writes into --out:

  mapping_qc_stats.tsv  every <dataset>.solo_qc.tsv, Dataset column prepended. The header
                        and row format are unchanged from the collectFile version, which
                        archive_run.sh, the resource notebook and triage all read.
                        Absent when no dataset was QC'd (e.g. --metaonly).
  skipped_samples.tsv   samples dropped for want of a species. Absent when none were.
  versions.yml          every task's versions.yml merged, one entry per process, sorted.
  datasets.tsv          per dataset: samples requested / in metadata / resolved / skipped /
                        mapped, runs, data types, species.
  samples.tsv           per resolved sample: species, chemistry, whitelist, runs, origins.
  runs.tsv              per run in links.tsv: sample, dataset, data type, species.
  run_chemistry.tsv     per run: RENAME10XRUN's chemistry inference and read lengths.

Standard library only, so it runs in any python3 image and standalone:
  collate_run_tables.py --datasets batch.csv --links links/ --sample-x-run sxr/ \
      --chemistry chem/ --soloqc qc/ --versions vers/ --samples samples.json \
      --skipped skipped.json --out tables/
"""
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

QC_HEADER = ("Dataset\tSample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\t"
             "Med_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\t"
             "exon_u\tfull_u+m\tfull_u")


def files_in(directory, pattern):
    d = Path(directory) if directory else None
    return sorted(d.glob(pattern)) if d and d.is_dir() else []


def read_tsv(path):
    with open(path, newline="") as fh:
        return [row for row in csv.reader(fh, delimiter="\t") if row and any(f.strip() for f in row)]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def load_json(path):
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return []
    return json.loads(Path(path).read_text())


def parse_versions(paths):
    """versions.yml is two levels deep: "PROCESS":\\n    tool: version. No YAML library needed."""
    merged = defaultdict(dict)
    for p in paths:
        process = None
        for line in Path(p).read_text().splitlines():
            if not line.strip():
                continue
            if not line.startswith((" ", "\t")):
                process = line.strip().rstrip(":").strip('"\'')
            elif process and ":" in line:
                tool, version = line.strip().split(":", 1)
                merged[process][tool.strip()] = version.strip()
    return merged


def write_versions(path, merged):
    with open(path, "w") as fh:
        for process in sorted(merged):
            fh.write(f'"{process}":\n')
            for tool in sorted(merged[process]):
                fh.write(f"    {tool}: {merged[process][tool]}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", required=True, help="the --datasets table (TSV: dataset_id, sample_id)")
    ap.add_argument("--links", help="directory of links.tsv files (run, species, urls, type, sample)")
    ap.add_argument("--sample-x-run", help="directory of <dataset>.sample_x_run.tsv files")
    ap.add_argument("--chemistry", help="directory of <run>.chemistry.json files")
    ap.add_argument("--soloqc", help="directory of <dataset>.solo_qc.tsv files")
    ap.add_argument("--versions", help="directory of versions*.yml files")
    ap.add_argument("--samples", help="JSON list of resolved sample metas")
    ap.add_argument("--skipped", help="JSON list of skipped sample metas")
    ap.add_argument("--extra-versions", help="a versions.yml to merge in as well (this task's own)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    # --datasets: what was asked for. An empty sample_id means "the whole series".
    requested = {}
    with open(a.datasets, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ds = (row.get("dataset_id") or "").strip()
            if ds:
                ids = [s.strip() for s in (row.get("sample_id") or "").split(",") if s.strip()]
                requested[ds] = ids

    # FETCH10XMETA: dataset -> sample -> runs, and per run its species and data type
    dataset_samples = defaultdict(set)
    sample_dataset, run_sample = {}, {}
    for p in files_in(a.sample_x_run, "*.sample_x_run.tsv"):
        ds = p.name[: -len(".sample_x_run.tsv")]
        for row in read_tsv(p):
            sample = row[0].strip()
            dataset_samples[ds].add(sample)
            sample_dataset.setdefault(sample, ds)
            for run in (row[1] if len(row) > 1 else "").split(","):
                if run.strip():
                    run_sample.setdefault(run.strip(), sample)
    run_species, run_type, sample_runs = {}, {}, defaultdict(set)
    for p in files_in(a.links, "*.tsv"):
        for row in read_tsv(p):
            if len(row) < 5:
                continue
            run, species, _urls, rtype, sample = (f.strip() for f in row[:5])
            run_species[run], run_type[run] = species, rtype
            run_sample.setdefault(run, sample)
            sample_runs[sample].add(run)

    # STARsolo QC, unchanged format
    qc_rows, mapped = [], defaultdict(set)
    for p in files_in(a.soloqc, "*.solo_qc.tsv"):
        ds = p.name[: -len(".solo_qc.tsv")]
        for row in read_tsv(p)[1:]:
            qc_rows.append([ds] + row)
            mapped[ds].add(row[0])
    if qc_rows:
        with open(out / "mapping_qc_stats.tsv", "w") as fh:
            fh.write(QC_HEADER + "\n")
            for row in qc_rows:
                fh.write("\t".join(row) + "\n")

    # Skipped samples, only when there are any
    skipped = load_json(a.skipped)
    if skipped:
        write_tsv(out / "skipped_samples.tsv", ["dataset_id", "sample_id", "specie"],
                  sorted([m.get("dataset_id", ""), m.get("id", ""), m.get("specie", "")] for m in skipped))
    skipped_by_ds = defaultdict(set)
    for m in skipped:
        skipped_by_ds[m.get("dataset_id", "")].add(m.get("id", ""))

    # Versions
    version_files = files_in(a.versions, "*.yml")
    if a.extra_versions:
        version_files.append(Path(a.extra_versions))
    write_versions(out / "versions.yml", parse_versions(version_files))

    # Samples that reached RENAME10XSAMPLE, with what was resolved for them
    samples = load_json(a.samples)
    resolved_by_ds = defaultdict(set)
    sample_rows = []
    for m in sorted(samples, key=lambda m: (m.get("dataset_id", ""), m.get("id", ""))):
        sid, ds = m.get("id", ""), m.get("dataset_id", "")
        resolved_by_ds[ds].add(sid)
        runs = sorted(sample_runs.get(sid, ()))
        sample_rows.append([
            sid, ds, m.get("specie", ""), m.get("chemistry", ""), m.get("wl", ""), len(runs),
            ",".join(sorted({run_type.get(r, "") for r in runs} - {""})),
            "yes" if sid in mapped.get(ds, ()) else "no",
        ])
    write_tsv(out / "samples.tsv",
              ["sample", "dataset", "specie", "chemistry", "wl", "n_runs", "origin", "mapped"], sample_rows)

    # Every run FETCH10XMETA listed, so a run that failed before RENAME10XRUN still maps to its sample
    write_tsv(out / "runs.tsv", ["run", "sample", "dataset", "type", "species"],
              sorted([r, run_sample.get(r, ""), sample_dataset.get(run_sample.get(r, ""), ""),
                      run_type.get(r, ""), run_species.get(r, "")] for r in run_type))

    # Per-run chemistry inference
    chem_rows = []
    for p in files_in(a.chemistry, "*.chemistry.json"):
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            print(f"warning: unreadable {p.name}", file=sys.stderr)
            continue
        run = d.get("run_id") or p.name[: -len(".chemistry.json")]
        sample = run_sample.get(run, "")
        lengths = d.get("read_lengths") or {}
        chem_rows.append([
            run, sample, sample_dataset.get(sample, ""),
            d.get("chemistry_id") or "", d.get("layout_confidence") or "", d.get("modality") or "",
            d.get("whitelist") or "", d.get("cb_len") or "", d.get("umi_len") or "",
            lengths.get("R1", ""), lengths.get("R2", ""), lengths.get("I1", ""),
            d.get("strand_hint") or "", d.get("ok"), len(d.get("warnings") or []),
        ])
    write_tsv(out / "run_chemistry.tsv",
              ["run", "sample", "dataset", "chemistry_id", "layout_confidence", "modality", "whitelist",
               "cb_len", "umi_len", "r1_len", "r2_len", "i1_len", "strand_hint", "ok", "n_warnings"],
              sorted(chem_rows, key=lambda r: (r[2], r[1], r[0])))

    # Per-dataset summary
    ds_rows = []
    for ds in sorted(set(requested) | set(dataset_samples)):
        in_meta = dataset_samples.get(ds, set())
        runs = {r for s in in_meta for r in sample_runs.get(s, ())}
        req = requested.get(ds)
        ds_rows.append([
            ds, len(req) if req else "all", len(in_meta), len(resolved_by_ds.get(ds, ())),
            len(skipped_by_ds.get(ds, ())), len(mapped.get(ds, ())), len(runs),
            ",".join(sorted({run_type.get(r, "") for r in runs} - {""})),
            ",".join(sorted({run_species.get(r, "") for r in runs} - {""})),
        ])
    write_tsv(out / "datasets.tsv",
              ["dataset", "requested", "in_metadata", "resolved", "skipped", "mapped", "runs",
               "data_types", "species"], ds_rows)

    print(f"{len(ds_rows)} datasets, {len(sample_rows)} samples, {len(chem_rows)} runs with chemistry, "
          f"{len(qc_rows)} QC rows, {len(skipped)} skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
