# Classifying and attributing failures

Distilled from `docs/agent_debug.md` §3, §4 and §6.

## §rules — the classifier

`bin/triage.py` holds the rule set; it is the single source of truth and this file does not
restate it. Read the `RULES` and `VERDICT` tables at the top of that file before arguing with a
category.

Two properties matter when you read its output:

* **First match wins, and order is load-bearing.** A block routinely matches several patterns.
  The order came from `scripts/triage3.py`; changing it changes the counts. Treat the buckets as
  a triage aid, not a taxonomy — if a number is going into a report, read the blocks behind it.
* **`UNCLASSIFIED` is where the interesting failures are.** In the August 2026 run, the two most
  valuable findings were sitting in a bucket of four: a silent single-mate SRA dump, and LSF
  memory kills wearing exit 130. `triage.py` prints every unclassified row with its work dir for
  exactly this reason. Read them; do not add a rule until you know what you are naming.

Adding a rule: put it above the more general pattern it must beat, give it a `VERDICT` entry in
the same commit, and re-run the batch-3 regression (`verify.md §regression`) to see what it
stole.

## §shape — two commands before anything else

```bash
# failures by process and exit code
awk -F'\t' 'NR>1{split($3,a," ");print a[1]"\t"$2}' data/tables/failedjobs6.tsv \
  | sort | uniq -c | sort -rn

# LSF kill reasons — these live only in the LSF report
grep -oE "TERM_[A-Z_]+" data/tables/failed6.log | sort | uniq -c
```

Then compare the shape against `known-issues.md §baseline`. A profile that differs sharply from
the August 2026 reference distribution is itself the finding.

## §exits — exit codes in this pipeline

| Process | Exit | Meaning |
|---|---|---|
| `RENAME10XRUN` | 2 | `TenxRunError` from `infer_10x_run_recommended.py`, **or** an argparse usage error. Both print to stderr; argparse dumps the whole option list, so look for the `error:` line at the end |
| `RENAME10XSAMPLE` | 2 | `RenameError` from `rename_fastqs_recommended.py` |
| `RENAME10XSAMPLE` | 1 | uncaught Python exception — batch 5's were `JSONDecodeError` on an empty `chemistry.json` |
| `WGET10X` | 8 | wget got an HTTP error response (403/404). **wget treats these as fatal** and ignores `--tries` unless `--retry-on-http-error` is set |
| `WGET10X` | 1 | wget usage error — in practice an empty URL field in the metadata |
| `SRA2FASTQ` | 1 | post-dump validation guard (see `verify.md §guards`) |
| `STARSOLO10X` | 1 | a check inside `starsolo_10x_auto.sh`, e.g. "No whitelist matched 200,000 random barcodes" |
| `STARSOLO10X` | 104, 139 | STAR's own fatal error. 104 has been the SJ-collapsed buffer; 139 is SIGSEGV |
| any | 130 | 128+2. In practice **LSF `TERM_MEMLIMIT`**, not a signal from the tool. Confirm with the `[lsf] TERM_*` line the collector appends, or `grep TERM_` in the raw log |
| any | 146, 140 | other LSF kills — `TERM_RUNLIMIT`, queue eviction |

**Never read `attempt` as diagnostic.** `nextflow.config` has

```groovy
errorStrategy = { task.exitStatus in ((130..145) + 104 + 175) && task.attempt < 3 || task.attempt == 1 ? 'retry' : 'ignore' }
```

`&&` binds tighter than `||`, so this parses as `(exitStatus in [...] && attempt < 3) || (attempt
== 1)` — **every** task retries once regardless of exit code. That is why 1795 of 1814 failures
in August 2026 sat at `attempt=2`. It is a property of the config, not a signal about the
failure. Still unfixed.

`errorStrategy 'ignore'` also means failures never reach the pipeline's exit status. A run can
finish `OK` in `.nextflow/history` with hundreds of failed tasks — batch 6 did.

## §attribution — run → sample → dataset

Tags are run accessions (`SRR…`) for `RENAME10XRUN` / `WGET10X` / `SRA2FASTQ`, and sample
accessions (`GSM…`) for `RENAME10XSAMPLE` / `STARSOLO10X`. `bin/triage.py` resolves both, but
the joins are worth knowing because they are where attribution silently fails:

```python
run2sample = {}                      # searchlist.tsv has NO header
for line in open("data/tables/searchlist.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 5: run2sample[p[0]] = p[4]

sample2dataset = {}                  # allhumandatasets FIRST: 4602 datasets vs datasets.tsv's 266
for fn in ("data/tables/allhumandatasets.tsv", "data/tables/datasets.tsv"):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        for s in r["sample_id"].split(","):
            sample2dataset.setdefault(s.strip(), r["dataset_id"])
```

Using `datasets.tsv` alone left 1336 of 1814 failures unattributed; with `allhumandatasets.tsv`
first, only 14 were. `searchlist.tsv` is per-batch and gets overwritten, so a high unattributed
count on a recent run usually means the searchlist is from a different batch — `triage.py` says
which files it used and how many rows it could not place.

**Always aggregate by dataset before reporting.** Failures are extremely concentrated: in August
2026 the top 5 datasets were 1495 of 1814 (82%), and GSE109816 alone was 48%. A flat list of
accessions hides that completely, and a fix aimed at the flat list aims at the wrong thing.

## §console — the live driver log

`logs/lsf/reprocessOutput<JOBID>.log` is the driver's stdout, ~1 MB per run. Line 2 gives the
run name. Failures appear as:

```
[ef/4eb2e2] NOTE: Process `…:SRA2FASTQ (SRR14567258)` terminated with an error exit status (1) -- Error is ignored
[50/635cdb] NOTE: Process `…:SRA2FASTQ (SRR14567256)` terminated with an error exit status (1) -- Execution is retried (1)
```

```bash
# per-process failure counts in a live run
grep -oE "terminated with an error exit status \([0-9]+\)" logs/lsf/reprocessOutput692614.log \
  | sort | uniq -c

# accession + exit code + work-dir hash, deduplicated
grep -oE "\[[0-9a-f]{2}/[0-9a-f]{6}\] NOTE: Process \`[^']*\([A-Z0-9]+\)\`.*status \([0-9]+\)" \
  logs/lsf/reprocessOutput692614.log | sort -u
```

`-- Error is ignored` means the task exhausted its retries and the pipeline moved on;
`-- Execution is retried (N)` means another attempt followed. The bracketed prefix resolves to
the work dir: `ls -d nf-work/ef/4eb2e2*`.

Use this only while a run is in flight, or when the cache is gone. Once the run has finished,
`bin/collect_run_logs.sh` gives you the same failures with exit codes, work dirs, retry counts
and recovery status already joined.

## §recovered — permanent versus self-healed

`failures<N>.tsv` carries a `recovered` column: did any attempt of that task name later reach
`COMPLETED` or `CACHED`? `triage.py` reports permanent and recovered separately and bases every
percentage on the permanent set.

This matters because the retry-everything-once bug means transient failures are always in the
raw failure list. Batch 3 had 111 failed tasks of which 19 self-healed; reporting 111 as the
failure count would have overstated the damage by 17%.
