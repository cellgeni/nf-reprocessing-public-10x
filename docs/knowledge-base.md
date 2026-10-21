---
title: Knowledge base
nav_order: 6
has_children: true
---

# Knowledge base

What public data does when you reprocess it at scale: the problems found in real submissions,
what each failure looks like, and the write-up of every past run. The rest of this site covers
how to run the pipeline; these pages cover what happens when you do.

| Page | Read it when |
|---|---|
| [Archive pathologies](archive-pathologies.md) | A dataset failed and you want to know whether it is a known problem, or you are about to conclude that a public submission is broken |
| [Reporting upstream](reporting-upstream.md) | You have confirmed a defect in a public submission and want to send it to GEO, SRA, ENA or ArrayExpress |
| [Failure modes](failure-modes.md) | You want to know what a given exit code means, what usually fails and why, or which problems are known and open |
| [tchem on our runs](tchem-benchmark.md) | You are deciding whether to use TusslChemistry (`tchem`), or its library-type and strand signals, beside our own chemistry inference. It records what tchem 0.1.0 did on 1,227 runs whose outcome we already knew |
| [Post-mortems](post-mortems.md) | You want to read the write-up of a past run — what failed, how much of it was the data and how much was us. Links to every published report |

The [10x chemistry reference](10x_chemistry_reference.md), which says what a set of FASTQs
actually is, sits under [Reference](reference/index.md).

## Why these pages exist

Public data is the hard part of this pipeline. The aligners are well understood and the failures
are rarely in them: they are in submissions that are mislabelled, internally inconsistent,
missing files, or perfectly correct but not what we assumed. The same small set of problems
recurs across batches, each one invisible until something several stages downstream fails
wearing an unrelated message.

Three habits are worth carrying into any use of these pages, because each one has been paid for:

- **Measure, do not infer — while you still can.** Several confident conclusions drawn from logs
  turned out to be wrong when checked against the FASTQs, because a failure routinely surfaces
  two stages after its cause. Work directories persist only until someone deletes them: batch 4
  was collected after its had gone, and all 47 of its failures now read `(no .command.log
  found)` — the counts survived in the trace, the reasons did not. Collect and archive a run
  while its work dirs still exist.
- **Aggregate by dataset before concluding anything.** Failures concentrate — one dataset was
  48% of a run's 1814 failures. A flat accession list hides that.
- **Separate "the data is broken" from "we did not check".** Most rejected volume in a typical
  run is correctly labelled data that we downloaded without reading its metadata. Confusing the
  two produces upstream reports that are wrong on their face and hides real pipeline bugs.

## For maintainers

Debugging a specific run — collecting its logs, classifying failures, attributing them to
datasets, verifying a chemistry claim against a work directory, writing the post-mortem — is
covered by the `reprocess-debug` skill in
[`.claude/skills/reprocess-debug/`](https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/.claude/skills/reprocess-debug),
which ships with the repository. [Developers](developers.md) describes the pipeline's
architecture and the traps in its channel wiring.

Complete records of past runs — the Nextflow trace, per-task tool stderr, execution reports and
the triage manifest — are archived on the Sanger farm at
`/nfs/cellgeni/reprocessing-runs/`, one directory per run. The readable write-up of each is
published as a browser-accessible page and indexed in [Post-mortems](post-mortems.md); the
farm archive is the durable copy, that index is where the links live.
