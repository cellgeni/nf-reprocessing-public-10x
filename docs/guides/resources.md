---
title: Resources and queues
parent: Guides
nav_order: 2
---

# Resources and queues
{: .no_toc }

1. TOC
{:toc}

## Labels

Every process carries a label for its size, as in nf-core. [`configs/base.config`](https://github.com/cellgeni/nf-reprocessing-public-10x/blob/main/configs/base.config)
says what each label means:

| Label | Resources | Processes |
|---|---|---|
| `process_single` | 1 cpu, 4 GB × attempt, 12 h | FETCH10XMETA, WGET10X, RENAME10XRUN, RENAME10XSAMPLE, REPROCESS10X_SRACHECK, REPROCESS10X_MAPPINGQC, REPROCESS10X_COLLATE, REPROCESS10X_RUNREPORT |
| `process_high` | 16 cpus | SRA2FASTQ, REPROCESS10X_SRAFIX, REPROCESS10X_BAM2FASTQ, STARsolo, Cell Ranger |
| `process_long` | 48 h | SRA2FASTQ, Cell Ranger |
| `process_network` | none: marks the tasks that need outbound HTTP | FETCH10XMETA, WGET10X |

## Per-process settings

What a label cannot say lives in `configs/<process>.config`:

| Process | Memory | Time | Other |
|---|---|---|---|
| FETCH10XMETA | 4 GB, +2 GB per retry | 12 h | at most 5 at once; up to 5 attempts |
| WGET10X | label | 12 h | at most 40 at once; up to 5 attempts; resumes in-job, verifies the archive and refetches it whole if the splice was bad |
| REPROCESS10X_SRACHECK | 16 GB × attempt | 1 h | |
| SRA2FASTQ | 16 GB, then 128 GB on any retry | 48 h | up to 5 attempts |
| REPROCESS10X_SRAFIX | 16 GB × attempt | 6 h × attempt | |
| REPROCESS10X_BAM2FASTQ | 4 GB × attempt | 12 h | |
| STARsolo | predicted per sample (below) | predicted | `--no-bam` |
| Cell Ranger | 128 GB, +64 GB per retry | 48 h | `--chemistry=auto --create-bam=false` |
| REPROCESS10X_MAPPINGQC | 2 GB, +4 GB per retry | 12 h | |
| REPROCESS10X_RUNREPORT | 8 GB × attempt | 12 h | |

The SRA2FASTQ rule comes from ~29,600 attempts: peak RSS was 4.1 GB median and 5.6 GB at p99, and
only two runs ever went past 8 GB (96 GB and 65 GB). More memory never fixed the tasks that failed
on the old 192-384 GB ladder; they failed deterministically.

The times are, for now, the run limit of the queue each process used on the farm before it had a
label (`normal` 12 h, `long` 48 h), so the farm routes exactly as it did. They are loose: the
longest successful runs over batches 20-26 took 11.9 h (RENAME10XRUN), 6.9 h (WGET10X), 4.0 h
(SRA2FASTQ) and 19.5 h (STARsolo).

## STARsolo

STARsolo's memory request and run limit are set for each sample from its input read count, which
is estimated before the task starts by reading the head of each R2 FASTQ
([`lib/StarsoloResources.groovy`](https://github.com/cellgeni/nf-reprocessing-public-10x/blob/main/lib/StarsoloResources.groovy)):

- **Memory**: (37.7 GB + 111.2 GB per 10⁹ reads) × 1.1, at least 48 GB, rounded up to 16 GB;
  ×1.5 on each retry, up to 1024 GB.
- **Run limit**: (0.13 h + 7.13 h per 10⁹ reads) × 2, in whole hours; doubled on each retry, up
  to 48 h.

The coefficients are q0.99 regressions fitted on batches 9-21 (11,998 samples) in
[`notebooks/resource_prediction.ipynb`](https://github.com/cellgeni/nf-reprocessing-public-10x/blob/main/notebooks/resource_prediction.ipynb).
On hold-out batches 22-24 they gave no memory or run-limit kills and sent about 95% of samples to
`normal`. Refit there and paste the printed constants into the class; do not hand-tune it. When
the estimate fails (no R2, an unreadable file), the old ladder applies: 128 GB × attempt, 48 h.

The run report compares the prediction with what each sample actually used — see
[Run report](../reference/run-report.md).

## Retries

Every task gets one retry, whatever its exit code. A third attempt follows only for the exit
codes known to be transient: 130-145 (signals, including the scheduler's memory kill and STAR's
segfault on a shared genome), 104 (STAR's own fatal error) and 175 (an LSF requeue). After that a
task is **ignored**, not fatal, so the rest of the run carries on; the run report lists every
sample lost that way. FETCH10XMETA, WGET10X and SRA2FASTQ have their own rule instead: up to 5
attempts whatever the code.

The first, unconditional retry is deliberate: transient failures are not identifiable from an
exit code. In one batch a STARsolo task died at exit 1 on a nonsense read length off a partly
staged FASTQ and succeeded on its second attempt. [Failure modes](../failure-modes.md) has the
codes and what usually causes them.

## Choosing the queue

The queue is the site profile's business, decided from the resolved memory and time; see
[How resources become queues](profiles.md#how-resources-become-queues). To run on another
cluster, write the same kind of rule for its queues.
