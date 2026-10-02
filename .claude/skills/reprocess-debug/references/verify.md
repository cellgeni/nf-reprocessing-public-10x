# Verifying a claim before you believe it

Distilled from `docs/agent_debug.md` §5, §8, §9 and §10.

**Work dirs persist until someone deletes them, and someone has.** `cleanup = false` in
`nextflow.config`, so Nextflow keeps them — but **the work dirs for batches 1-5 were deleted
(confirmed 2026-09-14)**, and batches 6-19 followed by 2026-10-02 (0 of their 818 failed-task
dirs remain, and batches 20-21's 134 went later the same day). Everything in this file works only on runs still on
scratch; for the rest there is nothing left to walk back to.

Check first, and do not promise a measurement you cannot take:

```bash
ls -d "$(awk -F'\t' 'NR==2{print $6}' data/tables/failures6.tsv)" 2>/dev/null \
  && echo "work dir present — measure it" || echo "gone — the archive is all there is"
```

Where they do survive this is the highest-leverage fact in the skill: **do not infer read
structure from log text when you can measure it.** Several confident-looking log-based
conclusions were wrong until checked. Where they do not, say that a claim rests on the archived
logs alone — batch 4 was collected after its work dirs went and all 47 of its failures carry
`(no .command.log found)`, so its counts are known and its reasons are not.

## §walkback — from a failed tag to the data

```bash
# work dir for a failed tag
awk -F'\t' '$2=="SRR11164904"{print $6}' data/tables/failures6.tsv

# staged inputs are symlinks into the UPSTREAM task's work dir — follow them to
# find that task's chemistry.json, .command.err, dump.log
for f in "$WD"/fastqs/*; do dirname "$(realpath "$f")"; done | sort -u
```

**A failure routinely surfaces two stages downstream of its cause, wearing an unrelated error
message.** `fastq-dump` dropping a mate appeared as an argparse usage error in `RENAME10XRUN`.
Truncated downloads appeared as gzip CRC errors. When a message looks like a tooling mistake
rather than a data problem, walk back up the staged symlinks before believing it.

## §whitelist — measuring the barcode match rate

Matching in both `infer_10x_run_recommended.py` and `starsolo_10x_auto.sh` is **exact set
membership on the first `cb_len` bases**. Consequences:

* A single `N` in the barcode is always a miss.
* The head of an Illumina FASTQ concentrates cycle-1 N-calls and the worst tile, so a head-only
  sample understates the match rate badly.

```python
import gzip, os
WLDIR = "/nfs/cellgeni/STAR/whitelists"
# v2 / 5' v1v2 = 737K-august-2016.txt (16bp) | v3+ = 3M-february-2018.txt (16bp)
# v1 = 737K-april-2014_rc.txt (14bp)         | arc = 737K-arc-v1.txt
# v4 3' = 3M-3pgex-may-2023.txt              | v4 5' = 3M-5pgex-jan-2023.txt
def probe(path, wl="737K-august-2016.txt", cb=16, start=0, n=100_000):
    """Match rate over a window, N-containing barcodes excluded from the denominator."""
    s = {l.strip() for l in open(os.path.join(WLDIR, wl))}
    tot = hasN = hit = 0
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            rec = i >> 2
            if rec < start: continue
            if rec >= start + n: break
            if i & 3 != 1: continue
            b = line.strip()[:cb]
            if len(b) < cb: continue
            tot += 1
            if "N" in b: hasN += 1
            elif b in s: hit += 1
    usable = tot - hasN
    return dict(n=tot, pct_N=100*hasN/tot if tot else 0,
                match=100*hit/usable if usable else 0)
```

Compare `start=0` against `start=1_000_000`. A large gap means a sampling artifact, not bad
data. A flat low rate in both, with low `pct_N`, means the data really is not 10x.

Read-geometry shorthand: `8+26+98` = I1 + (16 CB + 10 UMI, v2) + cDNA. `28` = 16+12, v3. `24` =
14+10, v1. Reverse-complement was checked on the August 2026 rejects and never helped; don't
spend time on it without a specific reason.

## §guards — current thresholds

`infer_10x_run_recommended.py` defaults relevant to triage:

| Flag | Default | Note |
|---|---|---|
| `--sample-records` | 200000 | reservoir sample size |
| `--sample-window` | 2000000 | records scanned; sample drawn uniformly. `0` = old prefix-only behaviour |
| `--min-whitelist-fraction` | 0.20 | denominator **excludes** N-containing barcodes |
| `--min-barcode-common-fraction` | 0.98 | constant-length bar; failing it falls through to the next row |
| `--min-barcode-usable-fraction` | 0.90 | fraction that must reach `cb_len` — the correctness bar |
| `--warn-barcode-umi-fraction` | 0.95 | below this, warn about reads STAR will drop for a short UMI |
| `--min-bio-usable-fraction` | 0.85 | was 0.98 |
| `--min-bio-read-length` | 40 | the real floor for a biological read, and what actually rejects a 20–25 bp ADT/HTO library |
| `--min-gex-r2-length` | 60 | **tie-breaker only** — applied when more than one candidate clears 40 bp. A single candidate below 60 bp is accepted, not rejected |

The whitelist summary table in the log has an `N-drop` column — sampled reads set aside for an N
in the barcode. A large `N-drop` with a low match rate is the head-sampling artifact, not bad
data.

The `--min-gex-r2-length` row is worth reading twice, because the table used to say it was "what
rejects ADT/HTO libraries" and that is wrong. In
`modules/cellgeni/rename10xrun/resources/usr/bin/infer_10x_run_recommended.py:956-963`:

```python
long_files = [i for i in remaining if i.common_length >= args.min_bio_read_length]   # 40
if allow_excluded:
    preferred = [i for i in long_files if i.common_length >= args.min_gex_r2_length]  # 60
    if len(preferred) == 1:                  # a lone sub-60bp candidate never gates
        long_files = preferred
bio = choose_one(long_files, "GEX biological R2 read")
```

A run with exactly one non-barcode read never reaches the 60 bp test: 20 bp and 25 bp reads are
rejected by the 40 bp floor, and a 42 bp read is accepted as GEX.

> **Do not try to separate feature-barcode libraries from GEX by read length.** Batch 9 measured
> both ends. GSE215253's two junk samples have a 42 bp biological read; `GSM6660161` (GSE216187,
> 56 bp), `GSM6668926` (GSE216329, 55 bp), `GSM6681083` (GSE216602, 55 bp) and `GSM6681046`
> (GSE216595, 55 bp) map at 78–94% with 1270–3165 median features. Raising the floor toward 60
> discards four good datasets to catch one bad pair. Geometry cannot do it either — GSE216914's
> CSP and VDJ libraries of one sample share an identical `10/10/26/90` layout and map at 0.9%
> and 93.9%. Use `§qc` instead.

## §qc — screening the outputs

Everything above this section is about samples that failed. This one is about samples that did
not, and it exists because batch 9 published 38 near-empty matrices with no task failing anywhere
and batch 10 published 75 — 4.8% and 10.0% of everything they aligned. `triage.py` cannot see
them; the only record is the QC table, and that is deleted with `results/` when the batch is
uploaded.

```bash
# de-duplicate on Sample first — see the trap below
awk -F'\t' 'NR>1 && !seen[$2]++ && $17+0 < 0.05 && $8+0 < 100 {print $1, $2, $14, $17, $8}' \
    results/batch10/mapping_qc_stats.tsv
```

Columns: `$1` Dataset, `$2` Sample, `$8` Med_nFeature, `$14` all_u+m, `$17` exon_u, `$19` full_u.

**Both terms are screening thresholds from two batches, not verdicts.** Treat a hit as "look at
this", not "this is junk".

* `exon_u < 0.05` does the work. Batch 9: 38 unique samples, maximum 0.0398, against a minimum of
  0.0597 across the other 746. Batch 10: 76.
* `Med_nFeature < 100` exists to spare single-nucleus data. On its own, `exon_u` called
  `GSM6729514` (GSE217892, batch 10) junk — `exon_u` 0.0425, but `all_u+m` 0.968, `full_u` 0.468
  and **583 median features**. snRNA-seq reads are intronic, so exonic is low while GeneFull is
  high. The second term drops that one sample and nothing else: batch 9 still flags all 38,
  batch 10 flags 75 rather than 76.
* **`full_u` is not the right second term**, tempting as it looks from that example. It would
  clear 7 of batch 9's 38 — GSE216999's cell-hashing libraries reach `full_u` 0.08-0.38 on 1-5
  median features. Tested and rejected; do not re-derive it.

**Corroborate with the test alignment, which is cheaper and earlier.** Before aligning, the
wrapper runs two 200,000-read test alignments to pick a strand and prints the result into the
STARsolo task's `.command.err`:

```
[WARN]  Low GeneFull mapping: forward=0%, reverse=0%
```

`max(forward, reverse) < 5%` flagged 31 of the 38 with zero false alarms. Note that 147 of batch
9's 917 completed tasks carry that warning at all, so the *warning* is not the signal — the
near-zero *value* is. This is also the number a future gate should use, since it is available
before the expensive alignment rather than after.

**Confirm what the library actually is from GEO, in one step.** A flagged sample is usually a
feature-barcode library, and the type is in the free-text title even though the structured
fields are useless:

```bash
grep -E "^!Sample_title|^!Sample_geo_accession" \
    results/batch9/metadata/GSE216999/GSE216999_family.soft | head
```

Look for `cell hashing`, `hashtag`, `HT`, `CSP`, `ADT`, `gRNA`, `enrichment PCR`,
`Custom library`, `VDJ`. Every one of batch 9's 38 is
`!Sample_library_strategy = RNA-Seq` and `!Sample_library_source = transcriptomic single cell`,
so the structured metadata will never filter them.

**What is not a hit.** Real data that merely maps low. GSE215908's 16 samples sit at `exon_u`
0.11–0.14 with 500–2300 median features and 2700–11500 cells — plausibly a xenograft or a
mislabelled organism, but usable data, and correctly untouched by the 0.05 threshold. Do not
widen the threshold to catch them; that is a different question, and
`docs/archive-pathologies.md` is where it is recorded.

### The other half: samples that never reached a task

The screen above only sees what was aligned. A sample dropped earlier leaves nothing at all — no
task, no work dir, no QC row. Compare what the batch asked for against what the run emitted:

```bash
python3 - <<'EOF'
import csv, os
B = "10"
# batch22+ live in batches_deduplicated/; the never-run batches/batch22-49 share the names
bt = next(p for p in (f"data/tables/batches_deduplicated/batch{B}.csv",
                      f"data/tables/batches/batch{B}.csv") if os.path.exists(p))
for r in csv.DictReader(open(bt), delimiter="\t"):
    ds  = r["dataset_id"]
    req = {s.strip() for s in r["sample_id"].split(",") if s.strip()}
    p   = f"results/batch{B}/metadata/{ds}/links.tsv"
    if not os.path.exists(p):
        print(f"{ds}: no links.tsv at all ({len(req)} samples)"); continue
    seen = {l.split("\t")[4].strip() for l in open(p) if len(l.split("\t")) >= 5}
    if req - seen:
        print(f"{ds}: {len(req & seen)}/{len(req)} emitted, missing {sorted(req - seen)}")
EOF
```

Batch 10 found GSE218936 emitting 2 of 8 requested samples with `FETCH10XMETA` exiting 0, because
all eight GSMs share two run sets of four. The survivor was then aligned against all four runs of
its group, so its matrix pools four experimental conditions (1G/µG × stim/unstim) under one
sample's name — at 96% mapping and 1865 median features, i.e. looking perfect. See
`issues/2026-09-19-batch10-metadata-drops-samples-silently.md`.

Two read-outs of that check are benign and should not be reported as losses:

* **`no links.tsv at all`** where the dataset's `FETCH10XMETA` task genuinely failed — it is
  already in the failure list. Batch 9's GSE215121 is this case.
* **A duplicate dataset row**, where the same samples are emitted under the pair's other
  accession. Check the sibling before calling anything lost.

Cell count is *not* a detector for the pooled case. Batch 10's two pooled matrices hold 32709 and
31193 cells against a batch median of 7471, but that is only about the 97th percentile and the
batch's largest legitimate sample has 70564.

Two operational traps for the first screen:

* **De-duplicate on `Sample` before counting anything.** The table carries one row per dataset
  unit, and duplicate dataset accessions mean a sample appears twice: batch 9, 917 rows for 784
  samples; batch 10, 973 for 748. Batch 10's 75 flagged samples first present as 97 rows.
* **Run this before the batch is cleaned up.** `results/<outdir>/` does not survive upload — in
  September 2026 only batches 9 and 10 still had one. `archive_run.sh` preserves
  `mapping_qc_stats.tsv` and the run's concatenated `links.tsv`, but only for runs archived
  after that was added; for batches 1–8 the output evidence is gone.

`SRA2FASTQ` fails rather than emitting a partial run. Its messages, all of which the classifier
recognises:

* `fastq-dump discarded a whole read … for having zero length in the archive` — the
  `Rejected N READS because READLEN < 1` case
* `fastq-dump produced N FASTQ file(s)` — fewer than 2 mates
* `mates of X disagree on record count` — partial dump
* `could not count records` / `not a whole number of FASTQ records` — ragged output

`WGET10X` rejects empty URLs, asserts exactly one file per task, and runs `pigz -t`/`gzip -t` on
the download. `ext.verify_downloads = false` disables the integrity test if transfer-queue time
becomes a problem.

`configs/starsolo10x.config` passes `--limitOutSJcollapsed 10000000` through `ext.args` after
`--`. The memory ladder does **not** help SJ-buffer failures — STAR fails deterministically
there, so every retry fails identically.

## §iterate — run the Python directly

Much faster than a pipeline run, and the only sane way to iterate on thresholds. Skip the
full-file record count, which takes minutes on multi-GB inputs:

```bash
python3 modules/cellgeni/rename10xrun/resources/usr/bin/infer_10x_run_recommended.py \
  --fastqs "$WD"/fastqs/* --run-id "$TAG" \
  --whitelist-dir /nfs/cellgeni/STAR/whitelists \
  --prefer-gex --ignore-mismatching-ids --no-check-counts --no-check-ids \
  --outdir . --json c.json --tsv c.tsv
```

**Stub-run a subworkflow** to validate channel wiring without touching data. Write the harness
outside the repo and cap resources — several processes request 16 CPUs:

```groovy
// local.config
process { executor='local'; cpus=1; memory='1 GB'; resourceLimits=[cpus:2, memory:'2 GB'] }
singularity.enabled = false
docker.enabled = false
```
```bash
nextflow -q run harness.nf -stub -c local.config -work-dir ./work
```

`REPROCESS10X_BAM2FASTQ`'s stub touches `<id>.fastq.gz` but declares `path("fastqs/*")`, so it
fails in stub mode — still broken, and it blocks stub-testing any BAM-origin path.

**Test bash inside a process by rendering it.** Extract the script block, substitute the
Nextflow placeholders, `bash -n`, then run against fixtures:

```python
src  = open("modules/cellgeni/sra2fastq/main.nf").read()
body = src.split('script:', 1)[1].split('"""', 2)[1]
body = (body.replace('${meta.id}', 'SRR9').replace('$task.cpus', '2')
            .replace('\\$', '$').split('cat <<-END_VERSIONS')[0])
open("body.sh", "w").write("#!/bin/bash -euo pipefail\n" + body)
```

**Fake a binary on PATH** to test the real rendered script end to end:

```groovy
process { beforeScript = 'export PATH=/path/to/fake/bin:$PATH' }
```

Used with a fake `wget` to exercise valid / truncated / empty-URL / 404 paths.

Bash pitfalls the generated scripts are prone to: a bare `wait` always returns 0 (so
backgrounded failures vanish), `$(cmd)` in an arithmetic context treats an empty result as 0,
and `set -u` kills the task on any typo'd variable.

## §cases — known-good regression accessions

Cheap to re-run after touching inference. Expected outcomes as of August 2026:

| Accession | What it is | Must |
|---|---|---|
| `SRR10742259` | v2 at 96.6%, R1 26 bp but only 68.1% constant | **pass** with a trimmed-barcode warning |
| `SRR11164904` | v2, 46% of barcodes contain N; 9.2% head → 78.8% windowed | **pass** |
| `SRR12273028` | 7.8% head → 80.3% at 1M reads | **pass** |
| `SRR6643897` | GSE109816; R1 = constant 11 bp `AACCAAGAGAT` + variable tail, 0.0% vs every whitelist both orientations | **stay rejected** — not 10x |
| `SRR11479076` | v3 barcode at 98%, R2 = 25 bp (ADT/HTO) | **stay rejected** |
| `GSM4227413` | 24 runs, R1 mixes 26/28 bp, all runs agree `gex_3pv2_or_5pv1v2` cb=16 umi=10 | **pass** with `--run-jsons` |
| `SRR10124100` / `SRR10124115` | `Rejected … READLEN < 1`, one mate only | `SRA2FASTQ` **must fail** |
| `SRR9895496` | ENA 403 that later served fine | transient — retry, don't blacklist |

Cost note: `--sample-window 2000000` takes ~41 s for a 3-file run vs ~23 s at 200000. Acceptable
on a 1-CPU task; relevant if you raise it much further.

## §regression — the classifier's own regression test

`scripts/triage3.py` is the frozen reference for batch 3. Any change to the rules in
`bin/triage.py` must leave its category counts untouched, or you have silently rewritten
history:

```bash
diff <(.claude/skills/reprocess-debug/bin/triage.py \
         --manifest reports/failures_2026-08-25_18-24-36.tsv \
       | sed -n '/process x exit/,/^$/p') \
     <(python3 scripts/triage3.py | sed -n '/process x exit/,/^$/p') \
  && echo identical
```

Second check, on a run with real Python tracebacks rather than a pre-extracted `first_error`:

```bash
.claude/skills/reprocess-debug/bin/triage.py \
  --failed-log data/tables/failed5.log --failedjobs data/tables/failedjobs5.tsv
```

Batch 5 should classify all 48 with 0 unattributed. If `UNCLASSIFIED` grows, read the blocks
before adding a rule.
