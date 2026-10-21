---
title: Run report
parent: Reference
nav_order: 4
---

# Run report
{: .no_toc }

`results/multiqc/multiqc_report.html` is built by REPROCESS10X_RUNREPORT, the last task of a run,
from the run's own Nextflow trace and the [run-level tables](outputs.md#run-level-tables).
Everything it plots is also in `results/multiqc/multiqc_data/` as tables.

1. TOC
{:toc}

## Where to start

**Samples without a mapped output.** Failed tasks are ignored rather than fatal, so a run can
finish cleanly and still have lost samples; this table is where they show up. It lists every
sample that did not end in a STARsolo output, with the furthest stage it reached and, for
failures, the process, exit class and work directory of the last failed attempt. It is the run's
failure manifest: [Failure modes](../failure-modes.md) says what each exit class usually means,
and [Archive pathologies](../archive-pathologies.md) whether a dataset is a known problem.

## Sections

| Section | Shows |
|---|---|
| Task outcomes | Every task per process: right first time, failed then recovered on a retry, or never succeeded; on a `-resume`, cached. Since every task retries once whatever its exit code, a second attempt is not itself a sign of trouble |
| Compute per process | Core-hours spent, the share spent on attempts that failed, failed attempts by exit class, and peak memory against the median request |
| Sample status | Per dataset: mapped, skipped for want of a species, failed at a stage (the first with a task that never succeeded), or not aligned (no aligner requested, or it never got that far) |
| Samples without a mapped output | As above |
| Datasets | From FETCH10XMETA: samples requested (all = the whole series), found in the metadata, resolved, skipped and mapped, with runs, data types and species |
| Chemistry calls | Runs per inferred chemistry, per dataset; `run_chemistry.tsv` has the per-run detail |
| Chemistry call confidence | How many calls are unique and how many `layout_only` — the read layout is fixed but one whitelist covers several kits (see the [10x chemistry reference](../10x_chemistry_reference.md)) |
| Samples | What each sample resolved to before alignment: species, Cell Ranger chemistry, STARsolo whitelist, runs merged and where they came from |
| STARsolo mapping | Per sample: cells, reads, fraction of reads in cells, good barcodes, median genes and UMIs, saturation, unique gene mapping, whitelist |
| Cells vs reads | One point per sample, to spot outliers at a glance |
| STARsolo resources: predicted vs observed | The reads estimated from the FASTQ heads, the memory request and run limit derived from them, and what each sample actually used. More than one attempt means the first request or limit was not enough — see [Resources and queues](../guides/resources.md#starsolo) |
| Software versions | Every tool's version, per process |

## Exit classes

Exit codes are grouped for the report: 130 is a scheduler memory kill on LSF, 139 a segfault
(usually STAR), 140 and 146 other LSF kills, and 104 STAR's own fatal error. Anything else is
shown as its code. On another scheduler the signal codes mean the same signals, but the kills
behind them may differ.
