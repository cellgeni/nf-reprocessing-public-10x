# Verifying a claim before you believe it

Distilled from `docs/agent_debug.md` §5, §8, §9 and §10.

**The work dirs persist.** `cleanup = false` in `nextflow.config`, and all 1814 work dirs from
the August 2026 run were still on disk months later. This is the highest-leverage fact in the
skill: **do not infer read structure from log text when you can measure it.** Several
confident-looking log-based conclusions were wrong until checked.

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
| `--min-gex-r2-length` | 60 | what rejects ADT/HTO libraries |

The whitelist summary table in the log has an `N-drop` column — sampled reads set aside for an N
in the barcode. A large `N-drop` with a low match rate is the head-sampling artifact, not bad
data.

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
