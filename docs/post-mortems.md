# Post-mortems

Every published failure post-mortem for this pipeline, with its link. **This file is the index
of record for those links** — other pages cite an individual report where it is the evidence for
something, but this is the only complete list, and a report missing from it is findable only by
listing every artifact the publishing account owns and guessing from titles.

Each report is published as a Claude Artifact: a self-contained HTML page hosted on claude.ai,
readable in a browser rather than as a file on the farm. The write-up is built from
`bin/triage.py --json` output, so its numbers are the trace's numbers and not retyped ones.

> **Who can open these.** An artifact link resolves for the account that published it and for
> anyone that account has shared it with — it is not public. The copy that needs no claude.ai
> account is the `postmortem.html` snapshot inside the run's archive directory on the farm,
> listed in the Archive column below. Both are kept deliberately: the artifact is the readable
> one, the archive is the durable one.

## Published

| Report | Covers | Published | Archive |
|---|---|---|---|
| [Batch 9 Post-Mortem](https://claude.ai/artifact/7FWxkLzbvRiCBHnxpXf39r) | Batch 9, run `gigantic_hypatia` — 108 failed tasks, 102 permanent, 55 samples lost of 839; two thirds of the failures are correct rejections, and the expensive findings are elsewhere — 9 good samples lost to a STAR segfault decided by a one-point margin in a test alignment, and **38 samples that never failed** published as near-empty matrices because nothing gates a sample on whether it mapped | 2026-09-18 | not yet archived |
| [Batch 8 Post-Mortem](https://claude.ai/artifact/EbcT6NhroVDPYwLWtypFYm) | Batch 8, run `cheeky_cuvier` — 53 failed tasks, 11 permanent, 6 samples and 2 dataset accessions lost; 10 of the 11 are one submission that registered each library's barcode and cDNA reads as separate SRA runs, and a fifth of the alignment compute went on duplicate dataset accessions | 2026-09-16 | `batch8/cheeky_cuvier/postmortem.html` |
| [Batch 7 Post-Mortem](https://claude.ai/code/artifact/42c9dcfc-ea39-4a24-a425-c58ef3cadfb9) | Batch 7, run `elegant_lamarr` — 10 failed tasks, 9 permanent, zero datasets lost; 8 trace to one missing value in the truncated-UMI chemistry guard, 1 is a genuinely near-empty source FASTQ | 2026-09-15 | not yet archived |
| [Batch 6 Post-Mortem](https://claude.ai/code/artifact/6f5e2eaf-b5e4-4476-b074-79683eec2333) | Batch 6, run `spontaneous_ampere` — 16 failed tasks, 7 permanent; two datasets emptied, and 98% of the wasted compute in two samples | 2026-09-14 | `batch6/spontaneous_ampere/postmortem.html` |
| [Batch 5 Reprocessing Post-Mortem](https://claude.ai/code/artifact/9d667def-0edb-4278-9a6d-50bc49c92931) | Batch 5, run `tender_brattain` — 48 failed tasks, 48 permanent | 2026-09-08 | `batch5/tender_brattain/postmortem.html` |
| [Run 337224 Post-Mortem](https://claude.ai/code/artifact/7d89d7e0-705a-4c5b-bdb3-621615a9504c) | LSF job 337224 / Nextflow run `evil_volhard`, the REQ-74217 follow-up. Source of the GSE247111 shared-BioSample library merge | 2026-09-01 | not archived — ran in a different working directory |
| [Batch 3 Failure Triage](https://claude.ai/code/artifact/ef2992f4-b281-4942-8f17-8ceb4c90fbb5) | Batch 3, run `big_keller` — 111 failed tasks, 92 permanent | 2026-09-01 | `batch3/big_keller/postmortem.html` |
| [Run 673887 Post-Mortem](https://claude.ai/code/artifact/547e62d8-45bb-4aa6-ad6a-0a95cd9d5f36) | LSF job 673887 — GSE188429, GSE189141, GSE183068, GSE188823 | 2026-08-25 | not archived — ran in a different working directory |
| [REQ-74217 Post-Mortem](https://claude.ai/code/artifact/bcb35389-362e-4e89-9309-597e073aa66b) | REQ-74217, 21 datasets, mixed public and internal. Source of the GSE247111 shared-BioSample finding and the ARC-v1 chemistry bug | 2026-08-21 | not archived — ran in a different working directory |
| [Reprocessing Run Post-Mortem](https://claude.ai/code/artifact/e4947587-de6b-4793-bd66-23da91af549a) | The August 2026 baseline — 1814 failed tasks, the reference distribution quoted throughout these pages. Not attributable to a named run; it predates the first batch by two days | 2026-08-17 | `legacy/august-2026-baseline` — evidence only, no HTML snapshot |

Archive paths are relative to `/nfs/cellgeni/reprocessing-runs/`. Each run directory also holds
the trace, per-task tool stderr and triage manifest the report was written from, described by its
own `manifest.txt` — which since 2026-09-14 carries the report's URL too, so the archive says
where its readable copy is rather than only that one exists.

Four of the ten have an HTML snapshot in the archive; the rest do not, either because the run was
never archived here or — for the August baseline — because `archive_run.sh` only recognises the
`batch<N>`/job-id source names. Where the Archive column says otherwise, the artifact is the only
copy outside scratch.

## Runs with no post-mortem

Recorded so that an absence reads as an absence rather than as a page someone forgot to link.

| Run | Why there is none |
|---|---|
| Batch 1, `peaceful_feynman` | No failure counts recorded for it at all |
| Batch 2, `distraught_mercator` | 294 failed tasks, never triaged into a report |
| Batch 4, `nauseous_spence` | 47 failed, 44 permanent, but collected after its work directories were deleted — every task reads `(no .command.log found)`, so the counts survive and the reasons do not. Nothing left to write a report from |

## Publishing a new one

The procedure — how the page is built, the house style that keeps the family recognisable, and
the theme bug that makes a report unreadable rather than merely ugly — is
[`references/report.md`](../.claude/skills/reprocess-debug/references/report.md) in the
`reprocess-debug` skill. Two steps matter for this file:

1. Publish the page as an artifact. Re-publishing the same source file redeploys to the **same**
   URL, so a correction never needs a new link; from a different conversation the artifact's
   `url` has to be passed explicitly or a duplicate is created instead.
2. **Add the row here**, in the same change. This is the step that gets skipped, because once the
   page is published the work feels finished.

The HTML source stays in `data/` (gitignored, on scratch) so a later session can edit and
redeploy rather than rebuild, and `bin/archive_run.sh` copies it into the run's archive
directory. Neither of those is the record of where the page is — this file is.
