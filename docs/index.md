# Documentation

Reference material and accumulated operational knowledge for
[nf-reprocessing-public-10x](https://github.com/cellgeni/nf-reprocessing-public-10x), a Nextflow
pipeline that fetches public 10x datasets from GEO, SRA, ENA and ArrayExpress and reprocesses
them uniformly with STARsolo or Cell Ranger.

For installation and parameters, see the [README](../README.md). These pages cover what the
README cannot: what public data does when you run it at scale.

## Pages

| Page | Read it when |
|---|---|
| [10x chemistry reference](10x_chemistry_reference.md) | You need to know what a set of FASTQs actually is — which whitelist implies which chemistry, what read geometry each layout has, and which layouts must be rejected from a gene-expression path rather than guessed at |
| [Archive pathologies](archive-pathologies.md) | A dataset failed and you want to know whether it is a known problem, or you are about to conclude that a public submission is broken |
| [Reporting upstream](reporting-upstream.md) | You have confirmed a defect in a public submission and want to send it to GEO, SRA, ENA or ArrayExpress |
| [Failure modes](failure-modes.md) | You want to know what a given exit code means, what usually fails and why, or which problems are known and open |

## Why these pages exist

Public data is the hard part of this pipeline. The aligners are well understood and the failures
are rarely in them: they are in submissions that are mislabelled, internally inconsistent,
missing files, or perfectly correct but not what we assumed. The same small set of problems
recurs across batches, each one invisible until something several stages downstream fails
wearing an unrelated message.

Three habits are worth carrying into any use of these pages, because each one has been paid for:

- **Measure, do not infer.** Work directories persist. Several confident conclusions drawn from
  logs turned out to be wrong when checked against the FASTQs, because a failure routinely
  surfaces two stages after its cause.
- **Aggregate by dataset before concluding anything.** Failures concentrate — one dataset was
  48% of a run's 1814 failures. A flat accession list hides that.
- **Separate "the data is broken" from "we did not check".** Most rejected volume in a typical
  run is correctly labelled data that we downloaded without reading its metadata. Confusing the
  two produces upstream reports that are wrong on their face and hides real pipeline bugs.

## For maintainers

Debugging a specific run — collecting its logs, classifying failures, attributing them to
datasets, verifying a chemistry claim against a work directory, writing the post-mortem — is
covered by the `reprocess-debug` skill in
[`.claude/skills/reprocess-debug/`](../.claude/skills/reprocess-debug/), which ships with the
repository. `CLAUDE.md` at the repository root describes the pipeline's architecture and the
traps in its channel wiring.

Complete records of past runs — the Nextflow trace, per-task tool stderr, execution reports and
the triage manifest — are archived on the Sanger farm at
`/nfs/cellgeni/reprocessing-runs/`, one directory per run.
