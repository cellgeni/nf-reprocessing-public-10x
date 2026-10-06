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
| [Batch 22 Post-Mortem](https://claude.ai/artifact/KTkk1FFghqYvVbSUvcinY4) | Batch 22, run `gloomy_church` (LSF 290127), the first batch from the deduplicated to-do list. 1,842 samples requested, **4 lost, all recoverable**; 73 failed tasks, **29 permanent, 19 of them harmless** (index-only runs, one 2,686-spot run). **Zero duplicate alignments**, against 9-23% in batches 8-21. The losses: 3 samples in SRA runs whose read layout changes partway through (GSE249159, GSE248489), refused at `SRA2FASTQ` although the data is complete, and 1 dropped at exit 0 by the shared-BioSample bug, its run pooled into its sibling (GSE254170, fifth recurrence). 3 samples published from partial input: ENA's FASTQs for two GSE249894 runs are truncated at source and match ENA's own MD5, and a 16-minute ENA outage left GSE252187's `GSM7996301` without a run. **190 of 1,838 published matrices are not GEX**, 116 of them V(D)J; the title scan flagged 44 healthy GEX samples (GSE249313, `scRNA-seq and TCR profiling …`) and missed 20. 35 STARsolo memory kills cost 2,686 core-hours. Hold list in `data/tables/hold_batch22.tsv` (193 hold, 4 review) | 2026-10-04 | `batch22/gloomy_church/postmortem.html` |
| [Batch 20-21 Post-Mortem](https://claude.ai/artifact/61wEcdbyagU3Tc8ztKZ7cF) | Batches 20 and 21, runs `sleepy_curie` and `sick_galileo`. Two concurrent runs, 1,517 samples requested, **239 lost, 230 of them recoverable**; 134 failed tasks, **121 permanent**. **`FETCH10XMETA` dropped 159 of batch 21's samples at exit 0** (GSE246613 147, GSE247111 12) through the shared-BioSample fallback, the fourth recurrence of batch 10's bug and the largest, and published 112 pooled matrices, most under a TCR or CellPlex GSM's name. **The STARsolo strand-test segfault cost 52 samples** (GSE245310): the test reads GeneFull, which on snRNA-seq calls 3′ data Reverse at margins of 9-23 points, so the drafted margin fix would not have helped. Split runs cost 15 more (GSE245998 split four ways and caught at `SRA2FASTQ`; GSE247205), now detectable from `sra.tsv` as runs of one SRX with equal spot counts. The unapplied truncated-UMI fix cost 4 (GSE247827). 175 published matrices are not GEX, 147 of them V(D)J. Batches 6-19's work dirs found deleted. Hold list in `data/tables/hold_batch20-21.tsv` (288 samples) | 2026-10-02 | `batch20/sleepy_curie/postmortem.html`, `batch21/sick_galileo/postmortem.html` |
| [Batch 18-19 Post-Mortem](https://claude.ai/artifact/2QVKvaZucsWTN1gNchQNF9) | Batches 18 and 19, runs `compassionate_shaw` and `focused_goldberg`. Two concurrent runs sharing several series under different accessions. 1,672 samples requested, 9 distinct lost; 51 failed tasks, **32 permanent, 27 of them harmless** (index-only runs and 0-spot SRA runs). The losses: batch 10's `FETCH10XMETA` shared-BioSample drop **recurred** (GSE241292: 5 lost, 5 published matrices pooling two libraries); **ENA's two-file FASTQ omitted the barcode read** that SRA holds (3 samples in GSE241998/9), and a retrospective scan showed batch 14's GSE228428 (5 samples) was the same defect, misfiled as a submitter one; and STARsolo discarded a finished 106k-cell matrix when Velocyto segfaulted on an 8.06 × 10⁹-read sample, on 6 attempts. **111 published matrices are not GEX**, 47 of them V(D)J libraries no screen catches. Hold list in `data/tables/hold_batch18-19.tsv` | 2026-09-30 | `batch18/compassionate_shaw/postmortem.html`, `batch19/focused_goldberg/postmortem.html` |
| [Batch 16-17 Post-Mortem](https://claude.ai/artifact/9D19BmmWAueLCUjsjv54UT) | Batches 16 and 17, runs `sad_brenner` and `cranky_cuvier`. Two concurrent runs, 1,557 samples requested, 35 lost, 255 failed tasks of which 244 are permanent. **239 of those 244 are correct rejections** (191 are one Smart-seq2 sample; 24 are index-reads-only deposits). The damage is outside the failure list. `FETCH10XMETA` dropped 15 of GSE234714's 20 samples at exit 0 and pooled four replicates into each of 5 published matrices. That is batch 10's defect, **now root-caused** to a BioSample fallback that ignores GEO's per-GSM SRX, with a diff. **193 published matrices are not GEX** (156 of batch 16's 737, the worst proportion yet), and the `exon_u` screen catches 28 of them. 92 of the 193 are lineage-barcode amplicons, a new class. 6 samples were published from partial input. Hold list in `data/tables/hold_batch16-17.tsv` | 2026-09-27 | `batch16/sad_brenner/postmortem.html`, `batch17/cranky_cuvier/postmortem.html` |
| [Batch 14-15 Post-Mortem](https://claude.ai/artifact/Uo2cXbCxuMb9hKE9ZdtupE) | Batches 14 and 15, runs `maniac_keller` and `infallible_dalembert` — two concurrent runs, 1,630 samples requested, 36 lost, and 44 failed tasks of which **27 are permanent**. Among the smallest failure lists recorded, and again not where the damage is: **164 published matrices (96 of 814 and 68 of 780) are not gene expression libraries**, 61 of them one dataset's `prot` arm that the title screen had no token for. 16 of the 27 permanent failures are three defects diagnosed in earlier batches and never fixed — the batch 9 STAR strand-test segfault (7 samples of measurably good data), the batch 7 truncated-UMI guard (a 96.8% v3 match rejected at `len=27`, fix drafted 2026-09-15 and still unapplied), and, new here, `FETCH10XMETA` destroying all 21 samples of GSE229166 for one unresolvable GSM | 2026-09-25 | `batch14/maniac_keller/postmortem.html`, `batch15/infallible_dalembert/postmortem.html` |
| [Batch 12-13 Post-Mortem](https://claude.ai/artifact/RHgyUCfAW6vPQvPFacxhn4) | Batch 12-13, run `boring_hawking` — 243 failed tasks, **65 permanent**, 27 samples correctly refused (24 of them hashtag libraries), and the batch's real cost elsewhere: **5 samples published from partial input**, `GSM7056649` built from 1 of its 17 runs and looking entirely healthy. 20 of the 65 are one missing value in `INDEX_LENGTHS`; 148 of 1441 aligned samples (10.3%) are not gene expression; batch 10's silent `FETCH10XMETA` drop did **not** recur | 2026-09-23 | `batch12-13/boring_hawking/postmortem.html` |
| [Batch 11 Post-Mortem](https://claude.ai/artifact/U5EBirz8SadJ6axeqzdUNb) | Batch 11, run `irreverent_elion` — the cleanest failure list of any batch (22 failed, **8 permanent**, no dataset lost) and the worst output quality measured: **144 of 844 aligned samples (17.1%) are not gene expression**. 97 are a bulk Ig-repertoire dataset admitted because chemistry inference has no upper bound on the barcode read length, 23 are V(D)J libraries that map at 73–90% and defeat every screen proposed so far, 24 are CITE/hashing. Batch 10's silent `FETCH10XMETA` drop did **not** recur | 2026-09-20 | `batch11/irreverent_elion` |
| [Batch 10 Post-Mortem](https://claude.ai/artifact/S85Ef5g7g6eiRyub9G19Gq) | Batch 10, run `trusting_carson` — the cleanest failure list in the series (16 failed, 13 permanent) and its worst outputs: **75 of 748 aligned samples** published as near-empty matrices, **6 samples dropped inside `FETCH10XMETA` at exit 0** with their reads pooled into two surviving matrices that look perfect, and 23.1% of alignment compute spent on duplicates | 2026-09-19 | `batch10/trusting_carson` |
| [Batch 9 Post-Mortem](https://claude.ai/artifact/7FWxkLzbvRiCBHnxpXf39r) | Batch 9, run `gigantic_hypatia` — 108 failed tasks, 102 permanent, 55 samples lost of 839; two thirds of the failures are correct rejections, and the expensive findings are elsewhere — 9 good samples lost to a STAR segfault decided by a one-point margin in a test alignment, and **38 samples that never failed** published as near-empty matrices because nothing gates a sample on whether it mapped | 2026-09-18 | `batch9/gigantic_hypatia` |
| [Batch 8 Post-Mortem](https://claude.ai/artifact/EbcT6NhroVDPYwLWtypFYm) | Batch 8, run `cheeky_cuvier` — 53 failed tasks, 11 permanent, 6 samples and 2 dataset accessions lost; 10 of the 11 are one submission that registered each library's barcode and cDNA reads as separate SRA runs, and a fifth of the alignment compute went on duplicate dataset accessions | 2026-09-16 | `batch8/cheeky_cuvier/postmortem.html` |
| [Batch 7 Post-Mortem](https://claude.ai/code/artifact/42c9dcfc-ea39-4a24-a425-c58ef3cadfb9) | Batch 7, run `elegant_lamarr` — 10 failed tasks, 9 permanent, zero datasets lost; 8 trace to one missing value in the truncated-UMI chemistry guard, 1 is a genuinely near-empty source FASTQ | 2026-09-15 | `batch7/elegant_lamarr/postmortem.html` |
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

Fourteen of the eighteen have an HTML snapshot in the archive (batches 9-11's directories hold a
`postmortem.html` although their rows name only the directory); the rest do not, either because
the run was never archived here or — for the August baseline — because `archive_run.sh` only
recognises the `batch<N>`/job-id source names. Where the Archive column says otherwise, the artifact is the only
copy outside scratch.

## Runs with no post-mortem

Recorded so that an absence reads as an absence rather than as a page someone forgot to link.

| Run | Why there is none |
|---|---|
| Batch 1, `peaceful_feynman` | Collected only on 2026-09-27, after its work dirs were deleted: 146 failed tasks, 137 permanent, but no failure reasons survive to write a report from |
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
