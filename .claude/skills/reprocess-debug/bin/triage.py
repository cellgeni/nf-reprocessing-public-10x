#!/usr/bin/env python3
"""Classify and attribute the failures of one reprocessing run.

    triage.py --manifest data/tables/failures6.tsv
    triage.py --failed-log data/tables/failed5.log --failedjobs data/tables/failedjobs5.tsv
    triage.py --manifest ... --json /tmp/triage6.json

The manifest form is preferred and is what bin/collect_run_logs.sh writes. The
failed-log form exists for runs collected before that script did, where the only
evidence is the block file.

Generalises scripts/triage3.py (batch 3, hardcoded). Rule order is first-match
wins and is taken from that script; changing the order changes the counts, so
treat the categories as a triage aid, not a taxonomy. If a count matters, read
the blocks.
"""
import argparse, collections, csv, json, os, re, sys

# --------------------------------------------------------------------- rules
# First match wins. Most specific first.
RULES = [
    ("meta-no-run-id",         r"No experiment or run ID found for"),
    ("read-counts-differ",     r"Read counts differ across selected run FASTQs"),
    ("chem-no-whitelist-hit",  r"No supported 10x run layout .*\(no whitelist hits\)"),
    ("chem-no-whitelist-rand", r"No whitelist matched .* random barcodes"),
    ("chem-ambiguous-layout",  r"No supported 10x run layout"),
    ("no-single-gex-r2",       r"Expected exactly one GEX biological R\d read"),
    ("bc-variable-length",     r"barcode read has variable lengths|has varying length"),
    ("r1-len-differs-runs",    r"read length differs across runs"),
    ("chem-group-mismatch",    r"Run-level chemistry metadata mismatch"),
    ("chem-json-unreadable",   r"JSONDecodeError|load_run_metadata"),
    ("truncated-gzip",         r"end-of-stream marker|CRC check failed|crc32 mismatch|not a valid (gzip|archive)"),
    ("corrupt-fastq",          r"FASTQ (header|plus line) at record .* does not start with"),
    ("sra-zero-length-read",   r"zero length in the archive|READLEN < 1"),
    ("sra-too-few-fastq",      r"expects 2-4 FASTQ files|fastq-dump produced \d+ FASTQ file"),
    ("sra-mate-mismatch",      r"disagree on record count"),
    ("sra-ragged-output",      r"could not count records|not a whole number of FASTQ records"),
    ("wget-empty-url",         r"no download URL|empty URL"),
    ("wget-exhausted",         r"did not complete after \d+ in-job attempts"),
    ("wget-http-error",        r"ERROR 40[0-9]|ERROR 50[0-9]"),
    ("star-sj-buffer",         r"buffer size for SJ output is too small"),
    ("star-segfault",          r"Segmentation fault"),
    ("lsf-memlimit",           r"TERM_MEMLIMIT"),
    ("lsf-runlimit",           r"TERM_RUNLIMIT"),
    ("star-fatal-other",       r"EXITING because of fatal error"),
    ("no-error-captured",      r"^\s*$|no \.command\.log found"),
]

# verdict, and the one-line reason it gets that verdict
VERDICT = {
    "meta-no-run-id":         ("investigate",       "GEO sample with no SRA experiment/run — metadata gap upstream"),
    "read-counts-differ":     ("pipeline-bug",      "ENA _1/_2/unsuffixed mix; wrong file set selected"),
    "chem-no-whitelist-hit":  ("investigate",       "no whitelist hit at all — check the window, not the head"),
    "chem-no-whitelist-rand": ("pipeline-bug",      "starsolo re-detect disagreeing with run-level inference"),
    "chem-ambiguous-layout":  ("too-strict",        "layout matched nothing the inference accepts"),
    "no-single-gex-r2":       ("correct-rejection", "feature-barcode (ADT/HTO) library, R2 too short for GEX"),
    "bc-variable-length":     ("too-strict",        "quality-trimmed barcode"),
    "r1-len-differs-runs":    ("pipeline-bug",      "missing chemistry JSON wiring"),
    "chem-group-mismatch":    ("pipeline-bug",      "runs of one sample disagree on chemistry"),
    "chem-json-unreadable":   ("pipeline-bug",      "empty or invalid chemistry.json reached RENAME10XSAMPLE"),
    "truncated-gzip":         ("infrastructure",    "download integrity"),
    "corrupt-fastq":          ("investigate",       "malformed record mid-file"),
    "sra-zero-length-read":   ("correct-rejection", "READLEN < 1"),
    "sra-too-few-fastq":      ("correct-rejection", "SRA2FASTQ guard fired"),
    "sra-mate-mismatch":      ("correct-rejection", "partial dump caught"),
    "sra-ragged-output":      ("correct-rejection", "ragged fastq-dump output caught"),
    "wget-empty-url":         ("pipeline-bug",      "empty URL field in the metadata"),
    "wget-exhausted":         ("infrastructure",    "ENA dropped transfer after in-job retries"),
    "wget-http-error":        ("infrastructure",    "HTTP 4xx/5xx — often transient, do not blacklist"),
    "star-sj-buffer":         ("pipeline-bug",      "raise --limitOutSJcollapsed; retries cannot help"),
    "star-segfault":          ("investigate",       "STAR SIGSEGV — check the reference and the read geometry"),
    "lsf-memlimit":           ("infrastructure",    "LSF killed the job, not a signal from the tool"),
    "lsf-runlimit":           ("infrastructure",    "LSF wall-clock kill"),
    "star-fatal-other":       ("investigate",       "STAR fatal error, read the block"),
    "no-error-captured":      ("evidence-lost",     "no stderr captured — nf-work cleaned?"),
    "UNCLASSIFIED":           ("investigate",       "read the block in full"),
}


def classify(text):
    return next((n for n, p in RULES if re.search(p, text)), "UNCLASSIFIED")


# --------------------------------------------------------------------- input
def load_manifest(path):
    rows = []
    for r in csv.DictReader(open(path, errors="replace"), delimiter="\t"):
        rows.append({
            "process": r.get("process", "?"),
            "tag": r.get("tag", "?"),
            "exit": r.get("exit", "?"),
            "recovered": r.get("recovered", "no"),
            "workdir": r.get("workdir", ""),
            "text": r.get("first_error") or "",
        })
    return rows


def load_blocks(path):
    """{tag: text} — task stderr only, LSF report stripped if still present."""
    blocks, cur, buf, inb = {}, None, [], False
    hdr = re.compile(r"^### (\S+) \(")
    for line in open(path, errors="replace"):
        m = hdr.match(line)
        if m:
            if cur:
                blocks[cur] = "\n".join(buf)
            cur, buf, inb = m.group(1), [], True
            continue
        if line.startswith("Sender: LSF System"):
            inb = False
        if inb:
            buf.append(line.rstrip())
    if cur:
        blocks[cur] = "\n".join(buf)
    return blocks


def load_failedjobs(path):
    out = {}
    for r in csv.DictReader(open(path, errors="replace"), delimiter="\t"):
        out[r["tag"]] = {
            "process": r["name"].split(" ")[0].split(":")[-1],
            "tag": r["tag"],
            "exit": r["exit"],
            "recovered": "no",
            "workdir": r.get("workdir", ""),
        }
    return out


# --------------------------------------------------------------- attribution
def base(tag):
    return tag.split(":")[0].replace("Loading ", "").strip()


def load_attribution(repo, searchlist, datasets):
    run2sample, run2type = {}, {}
    cands = [searchlist] if searchlist else [
        os.path.join(repo, "data/tables/searchlist.tsv"),
        os.path.join(repo, "data/searchlist.tsv"),
    ]
    used = []
    for fn in cands:
        if not os.path.exists(fn):
            continue
        used.append(fn)
        for line in open(fn, errors="replace"):          # NO header, 5 columns
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                run2sample.setdefault(p[0], p[4])
                run2type.setdefault(p[0], p[3])

    # allhumandatasets first: 4602 datasets vs datasets.tsv's 266. The other way
    # round left 1336 of 1814 failures unattributed.
    sample2dataset = {}
    for fn in (datasets or [os.path.join(repo, "data/tables/allhumandatasets.tsv"),
                            os.path.join(repo, "data/tables/datasets.tsv")]):
        if not os.path.exists(fn):
            continue
        used.append(fn)
        for r in csv.DictReader(open(fn, errors="replace"), delimiter="\t"):
            for s in (r.get("sample_id") or "").split(","):
                sample2dataset.setdefault(s.strip(), r["dataset_id"])
    return run2sample, run2type, sample2dataset, used


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", help="failures<N>.tsv from collect_run_logs.sh")
    ap.add_argument("--failed-log", help="failed<N>.log (block format)")
    ap.add_argument("--failedjobs", help="failedjobs<N>.tsv, required with --failed-log")
    ap.add_argument("--repo", default=os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))))
    ap.add_argument("--searchlist", help="override searchlist.tsv")
    ap.add_argument("--datasets", nargs="*", help="override dataset tables, in priority order")
    ap.add_argument("--top", type=int, default=20, help="datasets to list (default 20)")
    ap.add_argument("--json", help="write the classified rows and aggregates here")
    args = ap.parse_args()

    if args.manifest:
        rows = load_manifest(args.manifest)
        source = args.manifest
    elif args.failed_log and args.failedjobs:
        blocks = load_blocks(args.failed_log)
        jobs = load_failedjobs(args.failedjobs)
        rows = []
        for tag, j in jobs.items():
            rows.append(dict(j, text=blocks.get(tag, "")))
        missing = set(blocks) - set(jobs)
        if missing:
            print(f"note: {len(missing)} blocks with no failedjobs row", file=sys.stderr)
        noblock = set(jobs) - set(blocks)
        if noblock:
            print(f"note: {len(noblock)} failed tasks with no log block", file=sys.stderr)
        source = args.failed_log
    else:
        ap.error("give --manifest, or --failed-log with --failedjobs")

    if not rows:
        print("no failures")
        return 0

    run2sample, run2type, sample2dataset, used = load_attribution(
        args.repo, args.searchlist, args.datasets)

    def dataset_of(tag):
        t = base(tag)
        return sample2dataset.get(run2sample.get(t, t)) or sample2dataset.get(t) or "unresolved"

    for r in rows:
        r["cat"] = classify(r["text"])
        r["verdict"] = VERDICT.get(r["cat"], ("?", ""))[0]
        r["dataset"] = dataset_of(r["tag"])
        r["route"] = run2type.get(base(r["tag"]), "-")

    perm = [r for r in rows if r["recovered"] == "no"]
    recov = [r for r in rows if r["recovered"] != "no"]

    print(f"source               : {source}")
    print(f"attribution tables   : {', '.join(os.path.relpath(u, args.repo) for u in used) or 'none found'}")
    print(f"failed tasks         : {len(rows)}")
    print(f"  permanently failed : {len(perm)}")
    print(f"  recovered on retry : {len(recov)}")
    unres = sum(1 for r in perm if r["dataset"] == "unresolved")
    print(f"  unattributed       : {unres}"
          + ("   <- searchlist.tsv is from another batch, or these are sample-tagged"
             "\n                         tasks whose GSM is in no dataset table"
             if unres > len(perm) / 4 else ""))
    print()

    print("=== permanent failures: process x exit x category ===")
    print(f"{'process':<26}{'exit':<6}{'category':<24}n")
    for (p, e, c), n in collections.Counter(
            (r["process"], r["exit"], r["cat"]) for r in perm).most_common():
        print(f"{p:<26}{e:<6}{c:<24}{n}")
    print()

    print("=== permanent failures by category (with verdict) ===")
    cats = collections.Counter(r["cat"] for r in perm)
    vt = collections.Counter()
    for c, n in cats.most_common():
        v, note = VERDICT.get(c, ("?", ""))
        vt[v] += n
        print(f"  {c:<24}{n:>5}  {v:<18}{note}")
    print(f"  {'TOTAL':<24}{sum(cats.values()):>5}")
    assert sum(cats.values()) == len(perm), "categories must sum to the total"
    print()
    for v, n in vt.most_common():
        print(f"   {v:<18}{n:>5} ({100 * n / len(perm):.1f}%)")
    print()

    print(f"=== permanent failures by dataset (top {args.top}) ===")
    c = collections.Counter(r["dataset"] for r in perm)
    run = 0
    for i, (d, n) in enumerate(c.most_common()):
        run += n
        if i < args.top:
            print(f"  {d:<16}{n:>5}   cum {100 * run / len(perm):5.1f}%")
    print(f"  {len(c)} datasets in total")
    print()

    print("=== route x process (permanent) ===")
    for (rt, p), n in collections.Counter(
            (r["route"], r["process"]) for r in perm).most_common():
        print(f"  {rt:<8}{p:<26}{n}")
    print()

    print("=== distinct first errors, permanent, digits masked ===")
    msgs = collections.Counter(re.sub(r"\d[\d,._]*", "N", (r["text"] or "").strip()[:150])
                               for r in perm)
    for m, n in msgs.most_common(14):
        print(f"  {n:>4}  {m}")
    print()

    unc = [r for r in perm if r["cat"] == "UNCLASSIFIED"]
    if unc:
        print(f"=== UNCLASSIFIED ({len(unc)}) — read these blocks in full ===")
        for r in unc[:40]:
            print(f"  {r['tag']:<16}{r['process']:<26}{(r['text'] or '')[:110]}")
            print(f"      {r['workdir']}")
        if len(unc) > 40:
            print(f"  … and {len(unc) - 40} more")
        print()

    if recov:
        print("=== recovered on retry: what self-healed ===")
        for (p, rc), n in collections.Counter(
                (r["process"], r["cat"]) for r in recov).most_common():
            print(f"  {p:<26}{rc:<24}{n}")
        print()

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({
                "source": source,
                "total": len(rows),
                "permanent": len(perm),
                "recovered": len(recov),
                "by_category": dict(cats),
                "by_verdict": dict(vt),
                "by_dataset": dict(c),
                "by_process_exit_category": {
                    "|".join(k): v for k, v in collections.Counter(
                        (r["process"], r["exit"], r["cat"]) for r in perm).items()},
                "rows": rows,
            }, fh, indent=1)
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
