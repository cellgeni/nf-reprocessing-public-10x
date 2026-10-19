#!/usr/bin/env python3
"""One row per requested sample that did not end in a usable output, and why.

    sample_failures.py --batch 24 --run tender_torricelli \\
        --table data/tables/batches_deduplicated_sorted/batch24.csv --outdir results/batch24 \\
        --trace reports/execution_trace_2026-10-10_08-14-40.txt \\
        --manifest data/tables/failures24.tsv --out data/tables/sample_failures24.csv
    sample_failures.py --combine /nfs/cellgeni/reprocessing-runs

The first form writes sample_failures<N>.csv for one run. archive_run.sh calls it and
copies the result into the run's archive directory. The second form rebuilds
<root>/sample-failures.csv from every batch*/*/sample_failures*.csv under <root>, so
re-archiving a run never duplicates its rows. Two runs of one batch (a run and its
-resume) both stay in the combined file; that is history, not duplication. Filter to
the latest run of a batch for the current state.

Columns: see FIELDS. failedgroup is rerun / terminate / review. The reason vocabulary
and the category -> reason mapping live in triage.py (SAMPLE_REASON), next to the
rules, so the two cannot drift.

Work dirs are read while they exist (stderr text, chemistry.json), but no path is
written: nf-work is deleted after the debug session, so the evidence goes into
`comments` instead.

Velocyto is never checked. A sample whose Gene and GeneFull matrices are fine is
usable, whatever happened to Velocyto.
"""
import argparse, collections, csv, glob, json, os, re, sys

sys.dont_write_bytecode = True      # no __pycache__ inside the skill
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from triage import classify, load_blocks, load_manifest, SAMPLE_REASON  # noqa: E402

FIELDS = ["batch", "run", "dataset_id", "sample_id", "is10x", "modality", "failedgroup",
          "reason", "stage", "process", "exit", "category", "runs", "detected_by", "comments"]

STAGE = {"FETCH10XMETA": "metadata", "WGET10X": "download", "SRA2FASTQ": "download",
         "SRAFIX": "download", "BAM2FASTQ": "download", "RENAME10XRUN": "rename runs",
         "RENAME10XSAMPLE": "merge sample", "STARSOLO10X": "align", "CELLRANGER_COUNT": "align"}
RUN_LEVEL = ("WGET10X", "SRA2FASTQ", "SRAFIX", "BAM2FASTQ", "RENAME10XRUN")
SAMPLE_LEVEL = ("RENAME10XSAMPLE", "STARSOLO10X", "CELLRANGER_COUNT")
ALIGNERS = ("STARSOLO10X", "CELLRANGER_COUNT")
GROUP_RANK = {"terminate": 0, "rerun": 1, "review": 2}
REASON_RANK = {"pipeline bug": 0, "infrastructure": 1}

# Output screen (verify.md §qc). MEDIAN_GENES_MIN is a judgement call: batch 24's
# 28 non-GEX libraries had medians of 1-7 genes per cell, shallow real GEX sits in
# the hundreds.
MEDIAN_GENES_MIN = 20
ATAC_FRACTION_MIN = 0.05

# SKILL.md §4 title scan, token-delimited on purpose
TITLE_TOKENS = [
    ("vdj", re.compile(r"(^|[\s_-])(TCR|BCR|VDJ)([\s_-]|$)", re.I)),
    ("feature-barcode", re.compile(r"(^|[\s_-])(CITE|HTO|ADT|CSP|CMO|hashtag|hashing|gRNA|sgRNA)([\s_-]|$)", re.I)),
]


def one_line(s, n=300):
    s = re.sub(r"\s+", " ", (s or "")).strip()
    return s if len(s) <= n else s[:n - 1] + "…"


def short_process(process, name=""):
    """REPROCESS10X:STARSOLO10X_HUMAN:REPROCESS10X_STARSOLO10X -> STARSOLO10X."""
    p = process if process and process != "-" else (name or "").split(" (")[0]
    p = p.split(":")[-1]
    for k in STAGE:
        if p.endswith(k) or p == k:
            return k
    return p


def tag_id(tag):
    t = (tag or "").split(":")[0].strip()
    for pre in ("Converting .bam file to .fastq for ", "Loading "):
        if t.startswith(pre):
            t = t[len(pre):]
    return t.strip()


def read_table(path):
    """(dataset, sample) pairs requested. An empty sample_id asks for the whole series and
    is kept as (dataset, ""), for build() to expand into the samples metadata found."""
    out = []
    for r in csv.DictReader(open(path, errors="replace"), delimiter="\t"):
        samples = [s.strip() for s in (r.get("sample_id") or "").split(",") if s.strip()]
        out.extend((r["dataset_id"].strip(), s) for s in samples or [""])
    return list(dict.fromkeys(out))


def read_trace(path):
    """(process, id) -> attempts in trace order, each a dict."""
    tasks = collections.defaultdict(list)
    for r in csv.DictReader(open(path, errors="replace"), delimiter="\t"):
        p = short_process(r.get("process", ""), r.get("name", ""))
        if p not in STAGE:
            continue
        tasks[(p, tag_id(r.get("tag", "")))].append(r)
    return tasks


def outcome(attempts):
    st = [a.get("status", "") for a in attempts]
    if any(s in ("COMPLETED", "CACHED") for s in st):
        return "ok"
    if any(s == "FAILED" for s in st) and not any(s in ("ABORTED", "SUBMITTED", "RUNNING") for s in st[-1:]):
        return "failed"
    return "unfinished"


def read_file(path, tail=None):
    try:
        with open(path, errors="replace") as fh:
            txt = fh.read()
        return txt[-tail:] if tail else txt
    except OSError:
        return ""


def summary(path):
    d = {}
    for line in read_file(path).splitlines():
        k, _, v = line.partition(",")
        d[k] = v
    return d


def num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def atac_evidence(workdir, run):
    """Best ATAC whitelist fraction in a failed run's chemistry.json, if it survives."""
    try:
        js = json.load(open(os.path.join(workdir, f"{run}.chemistry.json")))
    except (OSError, ValueError):
        return None
    best = 0.0
    for st in js.get("stats") or []:
        for m in st.get("matches") or []:
            if m.get("modality") == "atac":
                best = max(best, num(m.get("fraction"), 0.0))
    return best


def soft_titles(outdir):
    titles = {}
    for f in glob.glob(os.path.join(outdir, "metadata", "*", "*_family.soft")):
        gsm = None
        for line in open(f, errors="replace"):
            if line.startswith("^SAMPLE = "):
                gsm = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_title") and gsm:
                titles[gsm] = line.split("=", 1)[1].strip()
    return titles


def title_modality(title):
    for mod, rx in TITLE_TOKENS:
        if rx.search(title or ""):
            return mod
    return None


def history_status(repo, run):
    st = None
    for line in read_file(os.path.join(repo, ".nextflow", "history")).splitlines():
        p = line.split("\t")
        if len(p) >= 4 and p[2] == run:
            st = p[3]
    return st


# --------------------------------------------------------------------- build
def build(args):
    requested = read_table(args.table)
    tasks = read_trace(args.trace)
    manifest = {}
    if args.manifest and os.path.exists(args.manifest):
        for r in load_manifest(args.manifest):
            manifest[(short_process(r["process"]), tag_id(r["tag"]))] = r
    # failed<N>.log keeps each failed task's whole stderr after nf-work is gone; the
    # manifest's first_error is one line and often not the line that names the cause
    blocks = {}
    if args.failed_log and os.path.exists(args.failed_log):
        for tag, text in load_blocks(args.failed_log).items():
            blocks.setdefault(tag_id(tag), text)

    links, have_links = collections.defaultdict(list), set()
    for f in glob.glob(os.path.join(args.outdir, "metadata", "*", "links.tsv")):
        ds = os.path.basename(os.path.dirname(f))
        have_links.add(ds)
        for line in open(f, errors="replace"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                links[(ds, p[4].strip())].append({"run": p[0], "species": p[1], "urls": p[2], "type": p[3]})

    # A whole-series row becomes every sample its metadata found; with no links.tsv it stays
    # (dataset, "") so the metadata branch below still reports the dataset as lost.
    expanded = []
    for ds, s in requested:
        found = sorted({k[1] for k in links if k[0] == ds}) if not s else []
        expanded.extend([(ds, x) for x in found] or [(ds, s)])
    requested = list(dict.fromkeys(expanded))

    skipped = {}
    sk = os.path.join(args.outdir, "skipped_samples.tsv")
    if os.path.exists(sk):
        for r in csv.DictReader(open(sk), delimiter="\t"):
            skipped[(r.get("dataset_id"), r.get("sample_id"))] = r.get("specie", "")

    qc = {}
    qf = os.path.join(args.outdir, "mapping_qc_stats.tsv")
    if os.path.exists(qf):
        for r in csv.DictReader(open(qf), delimiter="\t"):
            qc.setdefault((r.get("Dataset"), r.get("Sample")), r)

    titles = soft_titles(args.outdir)
    run_status = history_status(args.repo, args.run) or "?"

    rows = []

    def emit(ds, s, *, reason, group, stage, process="", exit_="", category="", runs="",
             detected_by, comments, modality="unknown", is10x="unknown"):
        hint = title_modality(titles.get(s, ""))
        if hint and modality in ("unknown", "feature-barcode") and group == "terminate":
            comments = f"{comments}; GEO title '{one_line(titles[s], 80)}' suggests {hint}"
            modality = hint
        rows.append({"batch": args.batch, "run": args.run, "dataset_id": ds, "sample_id": s,
                     "is10x": is10x, "modality": modality, "failedgroup": group, "reason": reason,
                     "stage": stage, "process": process, "exit": exit_, "category": category,
                     "runs": runs, "detected_by": detected_by, "comments": one_line(comments, 500)})

    def cause_for(proc, tid, attempts):
        """Classify one permanently failed task from its stderr."""
        m = manifest.get((proc, tid))
        last = attempts[-1]
        text = (m or {}).get("text", "")
        cat = classify(text) if text else "no-error-captured"
        wd = (m or {}).get("workdir") or last.get("workdir", "")
        how = "first error"
        full = (read_file(os.path.join(wd, ".command.err"), 20000) + "\n" +
                read_file(os.path.join(wd, ".command.log"), 20000)) if wd else ""
        if not full.strip():
            full = blocks.get(tid, "")
        # A full disk breaks whatever was writing, and the guard that trips next is not
        # the cause (batch 24: SRA2FASTQ's FASTQ-count check after pigz hit the quota).
        if cat != "disk-quota" and classify(full) == "disk-quota":
            cat, how = "disk-quota", "full stderr"
            text = next((l for l in full.splitlines() if re.search(r"Disk quota exceeded|No space left", l)), text)
        if cat in ("UNCLASSIFIED", "no-error-captured"):
            c2 = classify(full)
            if c2 not in ("UNCLASSIFIED", "no-error-captured"):
                cat, how = c2, "full stderr"
                hit = next((l for l in full.splitlines() if re.search(r"(?i)error|exceeded|fatal|abort", l)), "")
                text = hit or text
        reason, group, modality, is10x, caveat = SAMPLE_REASON[cat]
        if cat == "chem-no-whitelist-hit" and proc == "RENAME10XRUN" and wd:
            frac = atac_evidence(wd, tid)
            if frac is not None and frac >= ATAC_FRACTION_MIN:
                reason, group, modality, is10x = "unsupported modality", "terminate", "atac", "yes"
                caveat = f"ATAC whitelist matched {frac:.0%} of sampled barcodes"
        return {"proc": proc, "id": tid, "exit": last.get("exit", (m or {}).get("exit", "")),
                "cat": cat, "reason": reason, "group": group, "modality": modality,
                "is10x": is10x, "caveat": caveat, "text": text, "how": how}

    def check_output(ds, s):
        """("good"|"missing"|"bad", emit kwargs) for what is published under outdir."""
        base = os.path.join(args.outdir, "starsolo", ds, s)
        if not os.path.isdir(base):
            cands = glob.glob(os.path.join(args.outdir, "starsolo", "*", s))
            base = cands[0] if cands else base
        if not os.path.isdir(base):
            cr = glob.glob(os.path.join(args.outdir, "cellranger", "*", s, "**", "filtered_feature_bc_matrix*"),
                           recursive=True)
            return ("good", None) if cr else ("missing", None)
        o = os.path.join(base, "output")
        g = summary(os.path.join(o, "Gene", "Summary.csv"))
        uniq = num(g.get("Reads Mapped to Gene: Unique Gene"), 0.0)
        gene_mtx = glob.glob(os.path.join(o, "Gene", "filtered", "matrix.mtx*"))
        full_mtx = glob.glob(os.path.join(o, "GeneFull", "filtered", "matrix.mtx*"))
        if not gene_mtx or not full_mtx or uniq == 0:
            missing = [n for n, ok in (("Gene/filtered", gene_mtx), ("GeneFull/filtered", full_mtx)) if not ok]
            return "bad", dict(
                reason="pipeline bug", group="rerun", stage="output check", process="STARSOLO10X",
                exit_="0", detected_by="output check", modality="gex", is10x="yes",
                comments=("STARsolo exited 0 but the Gene matrix is unusable: "
                          f"Unique Gene {g.get('Reads Mapped to Gene: Unique Gene', 'n/a')}, "
                          f"genome unique {g.get('Reads Mapped to Genome: Unique', 'n/a')}, "
                          f"cells {g.get('Estimated Number of Cells', 'n/a')}"
                          + (f", missing {', '.join(missing)}" if missing else "")
                          + "; STAR does not check its temp-file writes, so a full disk produces this"))
        med = num(g.get("Median Gene per Cell"))
        q = qc.get((ds, s), {})
        exon_u, med_nf = num(q.get("exon_u")), num(q.get("Med_nFeature"))
        if (med is not None and med < MEDIAN_GENES_MIN) or \
           (exon_u is not None and med_nf is not None and exon_u < 0.05 and med_nf < 100):
            return "bad", dict(
                reason="unsupported modality", group="terminate", stage="output check",
                process="STARSOLO10X", exit_="0", detected_by="output check", is10x="yes",
                comments=(f"aligned, but median {g.get('Median Gene per Cell', 'n/a')} genes per cell, "
                          f"Unique Gene {g.get('Reads Mapped to Gene: Unique Gene', 'n/a')}, "
                          f"{g.get('Estimated Number of Cells', 'n/a')} cells: not a whole-transcriptome GEX "
                          "library (feature barcode, targeted panel or V(D)J)"))
        return "good", None

    for ds, s in requested:
        sample_runs = links.get((ds, s), [])
        aligner = [(p, tasks[(p, s)]) for p in ALIGNERS if (p, s) in tasks]
        aligned_ok = any(outcome(a) == "ok" for _, a in aligner)

        # --- 1. aligned in this run: check what it published
        if aligned_ok:
            if args.no_output_check:
                continue
            verdict, kw = check_output(ds, s)
            if verdict == "missing":
                emit(ds, s, reason="unknown", group="review", stage="output check",
                     process=aligner[0][0], exit_="0", detected_by="output check",
                     comments="aligner completed but no published output was found", modality="gex", is10x="yes")
            elif verdict == "bad":
                emit(ds, s, **kw)
            continue
        # Not aligned in this run, but a usable output is already published: an earlier
        # run of the session aligned it and this one re-ran it under a new task hash
        # (batch 24: meta.resources) and did not finish. Nothing is missing.
        if not args.no_output_check and check_output(ds, s)[0] == "good":
            continue

        # --- 2. metadata never produced links for the dataset
        if ds not in have_links:
            t = tasks.get(("FETCH10XMETA", ds))
            if t and outcome(t) == "failed":
                c = cause_for("FETCH10XMETA", ds, t)
                emit(ds, s, reason=c["reason"], group=c["group"], stage="metadata", process="FETCH10XMETA",
                     exit_=c["exit"], category=c["cat"], detected_by="task failure",
                     comments=f"whole dataset lost at metadata. {c['caveat']}. {c['how']}: {one_line(c['text'], 200)}")
            elif run_status == "ERR" or (t and outcome(t) == "unfinished"):
                emit(ds, s, reason="terminated run", group="rerun", stage="metadata", process="FETCH10XMETA",
                     detected_by="unfinished", comments=f"run {args.run} ended {run_status} before metadata was fetched")
            else:
                emit(ds, s, reason="unknown", group="review", stage="metadata", process="FETCH10XMETA",
                     detected_by="metadata", comments="no links.tsv for the dataset and no FETCH10XMETA failure recorded")
            continue

        # --- 3. dataset fetched, sample not in it
        if not sample_runs:
            t = tasks.get(("FETCH10XMETA", ds), [])
            err = read_file(os.path.join(t[-1].get("workdir", ""), ".command.err")) if t else ""
            if re.search(rf"no SRA/ENA metadata found for sample {re.escape(s)}\b", err):
                emit(ds, s, reason="no link available", group="terminate", stage="metadata",
                     process="FETCH10XMETA", detected_by="metadata",
                     comments="FETCH10XMETA found no SRA/ENA metadata for this sample and left it out")
            else:
                emit(ds, s, reason="pipeline bug", group="rerun", stage="metadata", process="FETCH10XMETA",
                     detected_by="metadata",
                     comments=("requested but not emitted by FETCH10XMETA, which exited 0; usually the "
                               "shared-BioSample drop that pools its reads into a sibling sample (#6)"
                               + ("" if err else "; FETCH10XMETA stderr no longer available")))
            continue

        # --- 4. no reference for the species
        if (ds, s) in skipped:
            sp = sorted({r["species"] for r in sample_runs if r["species"]})
            emit(ds, s, reason="no reference for species", group="terminate", stage="merge sample",
                 detected_by="species skip", is10x="yes",
                 comments=f"species resolved to UNKNOWN; run annotations: {', '.join(sp) or 'none'}")
            continue

        # --- 5. task failures on the way
        keys = [(p, r["run"]) for r in sample_runs for p in RUN_LEVEL if (p, r["run"]) in tasks]
        keys += [(p, s) for p in SAMPLE_LEVEL if (p, s) in tasks]
        causes, unfinished = [], []
        for k in keys:
            o = outcome(tasks[k])
            if o == "failed":
                causes.append(cause_for(k[0], k[1], tasks[k]))
            elif o == "unfinished":
                unfinished.append(k)
        renamed_ok = any(outcome(tasks[("RENAME10XRUN", r["run"])]) == "ok"
                         for r in sample_runs if ("RENAME10XRUN", r["run"]) in tasks)

        if causes:
            causes.sort(key=lambda c: (GROUP_RANK[c["group"]], REASON_RANK.get(c["reason"], 9)))
            c = causes[0]
            modality, is10x = c["modality"], c["is10x"]
            if renamed_ok and modality == "unknown":
                modality, is10x = "gex", "yes"
            failing_runs = sorted({x["id"] for x in causes if x["proc"] in RUN_LEVEL})
            others = collections.Counter(f"{x['proc']} {x['cat']}" for x in causes[1:])
            note = [c["caveat"], f"{c['how']}: {one_line(c['text'], 200)}"]
            if others:
                note.append("also: " + "; ".join(f"{k} x{n}" for k, n in others.items()))
            reason, group = c["reason"], c["group"]
            partial = c["proc"] in RUN_LEVEL and renamed_ok and len(failing_runs) < len(sample_runs)
            if partial:
                note.append(f"{len(failing_runs)} of {len(sample_runs)} runs failed, the rest passed chemistry")
                if reason in ("not 10x", "unsupported modality"):
                    # Leaving a non-GEX run out is correct, and the sample aligns from its
                    # GEX runs once the group is released at the end of the run. It has no
                    # output only because the run did not get that far, or not yet.
                    note.insert(0, f"its {c['cat']} run(s) will be left out")
                    reason, group = ("terminated run", "rerun") if run_status == "ERR" else ("unknown", "review")
                elif group == "terminate":
                    # The missing runs are real data that cannot be fetched or decoded, so the
                    # sample can only ever be published from part of its input (#12-13).
                    group = "review"
                    note.insert(0, "only partial input is available")
            emit(ds, s, reason=reason, group=group, stage=STAGE[c["proc"]], process=c["proc"],
                 exit_=c["exit"], category=c["cat"], runs=";".join(failing_runs), detected_by="task failure",
                 modality=modality, is10x=is10x, comments=". ".join(n for n in note if n))
            continue

        # --- 6. nothing failed, nothing finished
        modality, is10x = ("gex", "yes") if renamed_ok else ("unknown", "unknown")
        if unfinished or run_status == "ERR":
            k = unfinished[-1] if unfinished else None
            where = (f"{k[0]} {tasks[k][-1].get('status', '').lower()}" if k
                     else "the sample never got as far as a task")
            emit(ds, s, reason="terminated run", group="rerun",
                 stage=STAGE[k[0]] if k else "align", process=k[0] if k else "",
                 runs=k[1] if k and k[0] in RUN_LEVEL else "", detected_by="unfinished",
                 modality=modality, is10x=is10x,
                 comments=f"run {args.run} ended {run_status}; {where}")
        else:
            emit(ds, s, reason="unknown", group="review", stage="align", detected_by="unfinished",
                 modality=modality, is10x=is10x,
                 comments="no output and no failed task recorded for a run that finished OK")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    good = len(requested) - len(rows)
    print(f"{args.run}: {len(requested)} requested, {good} usable, {len(rows)} rows -> {args.out}", file=sys.stderr)
    for k, n in collections.Counter((r["failedgroup"], r["reason"]) for r in rows).most_common():
        print(f"  {n:>5}  {k[0]:<10} {k[1]}", file=sys.stderr)
    return 0


# ------------------------------------------------------------------- combine
def combine(root):
    def key(p):
        m = re.search(r"batch(\d+)", p)
        return (int(m.group(1)) if m else 10**9, p)
    files = sorted(glob.glob(os.path.join(root, "batch*", "*", "sample_failures*.csv")), key=key)
    out = os.path.join(root, "sample-failures.csv")
    tmp = out + ".tmp"
    n_all = 0
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for f in files:
            rd = csv.DictReader(open(f, newline="", errors="replace"))
            if rd.fieldnames != FIELDS:
                print(f"skip  {f}: header differs from the current FIELDS", file=sys.stderr)
                continue
            n = 0
            for r in rd:
                w.writerow(r)
                n += 1
            n_all += n
            print(f"{n:>6}  {os.path.relpath(f, root)}", file=sys.stderr)
    os.chmod(tmp, 0o664)        # anyone in the group archives runs, and rebuilds this
    os.replace(tmp, out)
    print(f"{n_all:>6}  rows -> {out}", file=sys.stderr)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combine", metavar="ROOT", help="rebuild ROOT/sample-failures.csv and exit")
    ap.add_argument("--batch")
    ap.add_argument("--run")
    ap.add_argument("--table", help="the batch's input table (TSV, comma-joined sample_id)")
    ap.add_argument("--outdir", help="results/batch<N>")
    ap.add_argument("--trace", help="reports/execution_trace_<stamp>.txt or runlogs<N>.tsv")
    ap.add_argument("--manifest", help="failures<N>.tsv from collect_run_logs.sh")
    ap.add_argument("--failed-log", help="failed<N>.log, full stderr per failed task (fallback once nf-work is gone)")
    ap.add_argument("--out", help="where to write sample_failures<N>.csv")
    ap.add_argument("--no-output-check", action="store_true",
                    help="trust a completed aligner task; use when outdir now holds a later run's outputs")
    ap.add_argument("--repo", default=os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 4)))
    args = ap.parse_args()
    if args.combine:
        return combine(args.combine)
    missing = [f for f in ("batch", "run", "table", "outdir", "trace", "out") if not getattr(args, f)]
    if missing:
        ap.error("missing --" + ", --".join(missing))
    return build(args)


if __name__ == "__main__":
    raise SystemExit(main())
