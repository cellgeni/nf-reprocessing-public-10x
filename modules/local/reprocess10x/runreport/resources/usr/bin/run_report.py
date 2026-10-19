#!/usr/bin/env python3
"""
Turn one run's trace, COLLATE tables and STARsolo outputs into MultiQC custom content.

Writes into --out (then run `multiqc --config <out>/multiqc_config.yaml <out>`):

  multiqc_config.yaml            title, run header, section order, table limits
  run_outcomes_mqc.json          per process: tasks right first time / recovered / never succeeded
  run_resources_mqc.json         per process: core-h, wasted core-h, exit classes, peak vs request
  sample_status_mqc.json         status of every sample the run touched, as counts per dataset
  sample_problems_mqc.json       every sample that did not end in a mapped STARsolo output, and why
  datasets_mqc.json              COLLATE's datasets.tsv
  chemistry_mqc.json             chemistry calls per dataset (RENAME10XRUN)
  layout_confidence_mqc.json     how sure each chemistry call was
  samples_mqc.json               COLLATE's samples.tsv (RENAME10XSAMPLE's resolved metadata)
  mapping_mqc.json               mapping_qc_stats.tsv + STARsolo Summary.csv in general stats
  cells_vs_reads_mqc.json        scatter
  resources_predicted_mqc.json   STARsolo predicted vs observed reads, memory and runtime
  software_mqc_versions.yml      COLLATE's merged versions.yml
  star/<sample>/Log.final.out    symlinks for MultiQC's own STAR module

JSON is valid YAML, so every file is plain json.dumps; standard library only. The task
logic (attempts grouped per task, recovered vs permanent, exit classes) is the same as
notebooks/resource_prediction.ipynb §2.

Standalone, e.g. against a live run:
  run_report.py --trace reports/execution_trace_<stamp>.txt --tables results/<run> \\
      --starsolo results/<run>/starsolo --out mqc/
"""
import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

_UNIT = {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1.0, "TB": 1e3}
_DUR = re.compile(r"([\d.]+)\s*(ms|d|h|m|s)")
_DUR_H = {"ms": 1 / 3.6e6, "s": 1 / 3600, "m": 1 / 60, "h": 1.0, "d": 24.0}

# Pipeline stages a sample passes through, in order, and the processes behind each.
STAGES = [
    ("download", ("WGET10X", "SRA2FASTQ", "SRAFIX", "BAM2FASTQ")),
    ("rename runs", ("RENAME10XRUN",)),
    ("merge sample", ("RENAME10XSAMPLE",)),
    ("align", ("STARSOLO10X", "CELLRANGER_COUNT")),
    ("qc", ("MAPPINGQC",)),
]
RUN_LEVEL = {"WGET10X", "SRA2FASTQ", "SRAFIX", "BAM2FASTQ", "RENAME10XRUN"}


def gb(s):
    if not s or s.strip() in ("-", ""):
        return None
    try:
        v, u = s.split()
        return float(v) * _UNIT[u]
    except (ValueError, KeyError):
        return None


def hours(s):
    if not s or s.strip() in ("-", ""):
        return None
    return sum(float(v) * _DUR_H[u] for v, u in _DUR.findall(s))


def exit_class(row):
    if row["status"] == "COMPLETED":
        return "ok"
    code = row.get("exit", "-")
    return {"130": "130 memlimit", "139": "139 segfault", "140": "140/146 LSF kill",
            "146": "140/146 LSF kill", "104": "104 STAR fatal"}.get(code, f"exit {code}")


def step_of(process):
    last = process.split(":")[-1]
    for _stage, procs in STAGES:
        for p in procs:
            if p in last:
                return p
    return last.replace("REPROCESS10X_", "")


TAG_PREFIXES = ("Converting .bam file to .fastq for ", "Converting sra file to .fastq for ",
                "Collecting QC stats for ", "Loading ")


def run_of(tag):
    """The accession in a task tag: `SRR1:file.fastq.gz` (WGET10X) and the prose tags of
    BAM2FASTQ and friends (`Converting .bam file to .fastq for SRR1`) both give SRR1."""
    t = tag.split(":", 1)[0].strip()
    for pre in TAG_PREFIXES:
        if t.startswith(pre):
            return t[len(pre):].strip()
    return t


def read_tsv_dicts(path):
    if not path or not Path(path).exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(out, name, obj):
    (out / name).write_text(json.dumps(obj, indent=1, default=str))


def section(id_, name, description, plot_type, data, **extra):
    d = {"id": id_, "section_name": name, "description": description, "plot_type": plot_type, "data": data}
    d.update(extra)
    return d


def load_trace(path):
    rows = read_tsv_dicts(path)
    tasks = defaultdict(list)
    for r in rows:
        if r.get("status") not in ("COMPLETED", "FAILED", "ABORTED", "CACHED"):
            continue
        r["step"] = step_of(r["process"])
        h = hours(r.get("realtime")) if r.get("realtime") not in (None, "-") else None
        r["hours"] = h if h is not None else (hours(r.get("duration")) or 0.0)
        try:
            r["core_h"] = float(r.get("cpus") or 1) * r["hours"]
        except ValueError:
            r["core_h"] = r["hours"]
        if r["status"] == "CACHED":
            r["core_h"] = 0.0  # paid for by the run that produced it (a -resume)
        r["req_gb"], r["peak_gb"] = gb(r.get("memory")), gb(r.get("peak_rss"))
        tasks[(r["process"], r["name"])].append(r)
    return rows, tasks


def task_outcome(attempts):
    if all(a["status"] == "CACHED" for a in attempts):
        return "cached"
    ok = any(a["status"] in ("COMPLETED", "CACHED") for a in attempts)
    failed = [a for a in attempts if a["status"] in ("FAILED", "ABORTED")]
    return "first try" if ok and not failed else ("recovered" if ok else "permanent")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace", required=True)
    ap.add_argument("--tables", required=True, help="directory holding COLLATE's tables")
    ap.add_argument("--starsolo", help="directory of STARsolo sample dirs, named by sample")
    ap.add_argument("--run-info", help="JSON: runName, sessionId, commitId, commandLine, start")
    ap.add_argument("--predicted", help="JSON list of [sample, reads, mem_gb, hours] from meta.resources")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tables = Path(a.tables)

    trace_rows, tasks = load_trace(a.trace)
    run_info = json.loads(Path(a.run_info).read_text()) if a.run_info and Path(a.run_info).exists() else {}

    # ---- run overview -------------------------------------------------------------------
    outcomes = {key: task_outcome(att) for key, att in tasks.items()}
    n_perm = sum(o == "permanent" for o in outcomes.values())
    total_core_h = sum(r["core_h"] for att in tasks.values() for r in att)
    wasted_core_h = sum(r["core_h"] for att in tasks.values() for r in att if r["status"] in ("FAILED", "ABORTED"))
    header = [
        {"Run": run_info.get("runName", "")},
        {"Session": run_info.get("sessionId", "")},
        {"Commit": run_info.get("commitId") or "unknown"},
        {"Started": str(run_info.get("start", ""))},
        {"Report written": datetime.now().strftime("%Y-%m-%d %H:%M")},
        {"Command": run_info.get("commandLine", "")},
        {"Tasks": f"{len(tasks):,} ({n_perm:,} never succeeded)"},
        {"Compute": f"{total_core_h:,.0f} core-h, {wasted_core_h:,.0f} on failed attempts "
                    f"({100 * wasted_core_h / total_core_h if total_core_h else 0:.1f}%)"},
    ]

    # ---- task outcomes and resources per process ---------------------------------------
    by_step = defaultdict(lambda: Counter())
    res = defaultdict(lambda: {"core_h": 0.0, "wasted": 0.0, "req": [], "peak": [], "exits": Counter()})
    for (process, _name), att in tasks.items():
        step = att[0]["step"]
        by_step[step][outcomes[(process, _name)]] += 1
        r = res[step]
        for x in att:
            r["core_h"] += x["core_h"]
            if x["status"] == "CACHED":
                continue
            if x["status"] == "COMPLETED":
                if x["peak_gb"] is not None:
                    r["peak"].append(x["peak_gb"])
                if x["req_gb"] is not None:
                    r["req"].append(x["req_gb"])
            else:
                r["wasted"] += x["core_h"]
                r["exits"][exit_class(x)] += 1
    write(out, "run_outcomes_mqc.json", section(
        "run_outcomes", "Task outcomes",
        "Every task of the run, per process: right first time, failed then recovered on a retry, "
        "or never succeeded; on a -resume, cached tasks were done by an earlier run. Every task retries once whatever its exit code (nextflow.config), so a "
        "second attempt is not itself a sign of trouble.",
        "bargraph", {s: dict(c) for s, c in sorted(by_step.items())},
        pconfig={"id": "run_outcomes_plot", "title": "Task outcomes per process", "ylab": "tasks",
                 "cpswitch_c_active": False},
        categories={"first try": {"color": "#1b9e77"}, "recovered": {"color": "#d95f02"},
                    "permanent": {"color": "#e7298a"}, "cached": {"color": "#999999"}}))

    def q(vals, p):
        if not vals:
            return None
        s = sorted(vals)
        return round(s[min(len(s) - 1, int(p * (len(s) - 1) + 0.5))], 1)

    res_rows = {}
    for step, r in sorted(res.items()):
        res_rows[step] = {
            "tasks": sum(by_step[step].values()),
            "core_h": round(r["core_h"], 1),
            "wasted_core_h": round(r["wasted"], 1),
            "wasted_pct": round(100 * r["wasted"] / r["core_h"], 1) if r["core_h"] else 0,
            "req_gb_median": q(r["req"], .5),
            "peak_gb_median": q(r["peak"], .5),
            "peak_gb_p99": q(r["peak"], .99),
            "failed_attempts": ", ".join(f"{k}: {v}" for k, v in r["exits"].most_common()) or "",
        }
    write(out, "run_resources_mqc.json", section(
        "run_resources", "Compute per process",
        "Core-hours spent, the share spent on attempts that failed, the failed attempts by exit class "
        "(130 is an LSF memory kill, 139 a STAR segfault, 140/146 other LSF kills), and peak memory "
        "against the median request on completed attempts.",
        "table", res_rows,
        pconfig={"id": "run_resources_table", "title": "Compute per process", "col1_header": "Process"},
        headers={
            "tasks": {"title": "Tasks", "format": "{:,.0f}"},
            "core_h": {"title": "Core-h", "format": "{:,.0f}"},
            "wasted_core_h": {"title": "Wasted core-h", "format": "{:,.0f}", "scale": "OrRd"},
            "wasted_pct": {"title": "Wasted %", "suffix": "%", "max": 100, "scale": "OrRd"},
            "req_gb_median": {"title": "Request GB", "description": "median memory request"},
            "peak_gb_median": {"title": "Peak GB", "description": "median peak RSS"},
            "peak_gb_p99": {"title": "Peak GB p99"},
            "failed_attempts": {"title": "Failed attempts by exit class"},
        }))

    # ---- sample status ------------------------------------------------------------------
    datasets_rows = read_tsv_dicts(tables / "datasets.tsv")
    samples_rows = read_tsv_dicts(tables / "samples.tsv")
    chem_rows = read_tsv_dicts(tables / "run_chemistry.tsv")
    qc_rows = read_tsv_dicts(tables / "mapping_qc_stats.tsv")
    skipped_rows = read_tsv_dicts(tables / "skipped_samples.tsv")

    runs_rows = read_tsv_dicts(tables / "runs.tsv")
    sample_dataset = {r["sample"]: r["dataset"] for r in samples_rows}
    run_sample = {}
    for r in runs_rows + chem_rows:
        if r.get("sample"):
            run_sample.setdefault(r["run"], r["sample"])
            sample_dataset.setdefault(r["sample"], r.get("dataset", ""))
    for r in skipped_rows:
        sample_dataset.setdefault(r["sample_id"], r["dataset_id"])
    mapped = {r["Sample"] for r in qc_rows}
    skipped = {r["sample_id"] for r in skipped_rows}

    # task outcomes per sample and stage; run-level tags resolve through run_sample
    sample_stage = defaultdict(lambda: defaultdict(list))
    for (process, name), att in tasks.items():
        step = att[0]["step"]
        stage = next((s for s, procs in STAGES if step in procs), None)
        if stage is None or stage == "qc":
            continue
        tag = att[0].get("tag", "")
        sample = run_sample.get(run_of(tag), run_of(tag)) if step in RUN_LEVEL else tag
        sample_stage[sample][stage].append(att)

    status_rows, problems = {}, {}
    for sample in sorted(set(sample_dataset) | set(sample_stage)):
        ds = sample_dataset.get(sample, "")
        stages = sample_stage.get(sample, {})
        reached, failure = "", None
        for stage, _procs in STAGES:
            if stage not in stages:
                continue
            outs = [task_outcome(att) for att in stages[stage]]
            if any(o != "permanent" for o in outs):
                reached = stage
            perm = [att for att in stages[stage] if task_outcome(att) == "permanent"]
            if perm and failure is None:
                last = perm[0][-1]
                failure = (stage, last["step"], exit_class(last), last.get("workdir", ""))
        if sample in mapped:
            status = "mapped"
        elif sample in skipped:
            status = "skipped (species)"
        elif failure:
            status = f"failed: {failure[0]}"
        else:
            status = "not aligned"
        status_rows.setdefault(ds or "unknown", Counter())[status] += 1
        if status != "mapped":
            problems[sample] = {
                "dataset": ds, "status": status, "reached": reached or "-",
                "process": failure[1] if failure else "", "exit": failure[2] if failure else "",
                "workdir": failure[3] if failure else "",
            }
    write(out, "sample_status_mqc.json", section(
        "sample_status", "Sample status",
        "What became of every sample per dataset: mapped by STARsolo, skipped for want of a "
        "species, failed at a stage (the first stage with a task that never succeeded), or not "
        "aligned (no aligner requested, or it never got that far).",
        "bargraph", {ds: dict(c) for ds, c in sorted(status_rows.items())},
        pconfig={"id": "sample_status_plot", "title": "Sample status per dataset", "ylab": "samples"}))
    write(out, "sample_problems_mqc.json", section(
        "sample_problems", "Samples without a mapped output",
        f"The {len(problems):,} samples that did not end in a STARsolo output, with the furthest "
        "stage they reached and, for failures, the process, exit class and work dir of the last "
        "failed attempt. This is the run's failure manifest; triage it with the reprocess-debug skill.",
        "table", problems,
        pconfig={"id": "sample_problems_table", "title": "Samples without a mapped output",
                 "col1_header": "Sample"}))

    # ---- datasets, chemistry, samples ---------------------------------------------------
    if datasets_rows:
        write(out, "datasets_mqc.json", section(
            "datasets", "Datasets (FETCH10XMETA)",
            "Per dataset: samples requested (all = the whole series), found in the archive metadata, "
            "resolved to a species, skipped and mapped, plus runs, data types and species.",
            "table", {r["dataset"]: {k: v for k, v in r.items() if k != "dataset"} for r in datasets_rows},
            pconfig={"id": "datasets_table", "title": "Datasets", "col1_header": "Dataset"}))
    if chem_rows:
        chem = defaultdict(Counter)
        conf = defaultdict(Counter)
        for r in chem_rows:
            chem[r.get("dataset") or "unknown"][r.get("chemistry_id") or "none"] += 1
            conf[r.get("dataset") or "unknown"][r.get("layout_confidence") or "none"] += 1
        write(out, "chemistry_mqc.json", section(
            "chemistry", "Chemistry calls (RENAME10XRUN)",
            "Runs per inferred chemistry, per dataset. Full per-run detail (read lengths, whitelist, "
            "strand, warnings) is in run_chemistry.tsv.",
            "bargraph", {k: dict(v) for k, v in sorted(chem.items())},
            pconfig={"id": "chemistry_plot", "title": "Chemistry per dataset", "ylab": "runs"}))
        write(out, "layout_confidence_mqc.json", section(
            "layout_confidence", "Chemistry call confidence",
            "layout_only calls fix the read layout but not a unique chemistry; see "
            "docs/10x_chemistry_reference.md.",
            "bargraph", {k: dict(v) for k, v in sorted(conf.items())},
            pconfig={"id": "layout_confidence_plot", "title": "Chemistry call confidence", "ylab": "runs"}))
    if samples_rows:
        write(out, "samples_mqc.json", section(
            "samples", "Samples (RENAME10XSAMPLE)",
            "What each sample resolved to before alignment: species, Cell Ranger chemistry, STARsolo "
            "whitelist, runs merged and where they came from.",
            "table", {r["sample"]: {k: v for k, v in r.items() if k != "sample"} for r in samples_rows},
            pconfig={"id": "samples_table", "title": "Samples", "col1_header": "Sample"}))

    # ---- mapping ------------------------------------------------------------------------
    star_dir = out / "star"
    summary = {}
    if a.starsolo and Path(a.starsolo).is_dir():
        for d in sorted(Path(a.starsolo).iterdir()):
            log = d / "Log.final.out"
            if log.exists():
                (star_dir / d.name).mkdir(parents=True, exist_ok=True)
                link = star_dir / d.name / "Log.final.out"
                if not link.exists():
                    os.symlink(os.path.realpath(log), link)
            s = d / "output" / "Gene" / "Summary.csv"
            if s.exists():
                kv = dict(line.split(",", 1) for line in s.read_text().splitlines() if "," in line)
                summary[d.name] = kv
    gen = {}
    for r in qc_rows:
        s = r["Sample"]
        gen[s] = {"cells": r.get("Cells"), "reads_m": round(float(r.get("Rd_all") or 0) / 1e6, 1),
                  "frc_in_cells": r.get("Frc_in_cells"), "good_bc": r.get("Good_BC"),
                  "med_nfeature": r.get("Med_nFeature"), "wl": r.get("WL")}
    for s, kv in summary.items():
        g = gen.setdefault(s, {})
        for key, col in (("Sequencing Saturation", "saturation"), ("Median UMI per Cell", "median_umi"),
                         ("Reads Mapped to Gene: Unique Gene", "gene_unique")):
            if key in kv:
                g[col] = kv[key].strip()
    if gen:
        write(out, "mapping_mqc.json", section(
            "mapping", "STARsolo mapping", "From mapping_qc_stats.tsv and STARsolo's Gene/Summary.csv.",
            "generalstats", gen,
            pconfig=[
                {"cells": {"title": "Cells", "format": "{:,.0f}"}},
                {"reads_m": {"title": "M reads", "description": "STAR input reads, millions"}},
                {"frc_in_cells": {"title": "Frac in cells", "max": 1, "min": 0}},
                {"good_bc": {"title": "Good BC", "max": 1, "min": 0}},
                {"med_nfeature": {"title": "Median genes", "format": "{:,.0f}"}},
                {"saturation": {"title": "Saturation", "max": 1, "min": 0}},
                {"median_umi": {"title": "Median UMI", "format": "{:,.0f}"}},
                {"gene_unique": {"title": "Unique gene", "max": 1, "min": 0,
                                 "description": "fraction of reads mapped uniquely to a gene"}},
                {"wl": {"title": "WL"}},
            ]))
        scatter = {s: {"x": float(r.get("Rd_all") or 0) / 1e6, "y": float(r.get("Cells") or 0)}
                   for s, r in ((r["Sample"], r) for r in qc_rows)}
        write(out, "cells_vs_reads_mqc.json", section(
            "cells_vs_reads", "Cells vs reads", "One point per sample.", "scatter", scatter,
            pconfig={"id": "cells_vs_reads_plot", "title": "Cells vs input reads",
                     "xlab": "input reads (M)", "ylab": "cells"}))

    # ---- predicted vs observed STARsolo resources --------------------------------------
    predicted = json.loads(Path(a.predicted).read_text()) if a.predicted and Path(a.predicted).exists() else []
    observed_reads = {r["Sample"]: float(r.get("Rd_all") or 0) for r in qc_rows}
    star_done = {}
    for (process, _n), att in tasks.items():
        if att[0]["step"] == "STARSOLO10X":
            done = [x for x in att if x["status"] in ("COMPLETED", "CACHED")]
            star_done[att[0].get("tag", "")] = (done[-1] if done else None, len(att))
    pred_rows = {}
    for row in predicted:
        if not row or row[1] is None:
            continue
        sample, reads, mem_gb, limit_h = row[:4]
        obs, n_att = star_done.get(sample, (None, 0))
        true_reads = observed_reads.get(sample)
        pred_rows[sample] = {
            "pred_reads_m": round(reads / 1e6, 1),
            "obs_reads_m": round(true_reads / 1e6, 1) if true_reads else None,
            "reads_err_pct": round(100 * (reads / true_reads - 1), 1) if true_reads else None,
            "mem_req_gb": mem_gb, "peak_gb": round(obs["peak_gb"], 1) if obs and obs["peak_gb"] else None,
            "limit_h": limit_h, "realtime_h": round(obs["hours"], 2) if obs else None,
            "attempts": n_att,
        }
    if pred_rows:
        write(out, "resources_predicted_mqc.json", section(
            "resources_predicted", "STARsolo resources: predicted vs observed",
            "Reads estimated from the FASTQ heads and the memory request and run limit derived from them "
            "(lib/StarsoloResources.groovy), against STAR's input reads, the peak RSS and the runtime. "
            "More than one attempt means the first request or limit was not enough.",
            "table", pred_rows,
            pconfig={"id": "resources_predicted_table", "title": "STARsolo resources", "col1_header": "Sample"}))

    # ---- versions and config ------------------------------------------------------------
    if (tables / "versions.yml").exists():
        (out / "software_mqc_versions.yml").write_text((tables / "versions.yml").read_text())

    order = ["run_outcomes", "run_resources", "sample_status", "sample_problems", "datasets",
             "chemistry", "layout_confidence", "samples", "star", "cells_vs_reads", "resources_predicted"]
    config = {
        "title": f"Reprocessing run {run_info.get('runName', '')}".strip(),
        "subtitle": "nf-reprocessing-public-10x",
        "report_header_info": header,
        "report_section_order": {sid: {"order": 1000 - 10 * i} for i, sid in enumerate(order)},
        "output_fn_name": "multiqc_report.html",
        "data_dir_name": "multiqc_data",
        "max_table_rows": 5000,
        "plots_force_flat": False,
        "plots_flat_numseries": 2000,
        "sp": {"star": {"fn": "*Log.final.out"}},
        "show_analysis_paths": False,
        "no_version_check": True,   # execution nodes may have no route out; it only costs a timeout
    }
    (out / "multiqc_config.yaml").write_text(json.dumps(config, indent=1, default=str))

    print(f"{len(tasks):,} tasks, {len(status_rows)} datasets, {len(problems):,} samples without output, "
          f"{len(gen):,} mapped samples, {len(pred_rows):,} with predictions", file=sys.stderr)


if __name__ == "__main__":
    main()
