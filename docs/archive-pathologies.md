# Archive pathologies

A catalogue of defects and traps found in public 10x submissions while reprocessing them, one
row per accession, each backed by evidence still on disk.

It exists because these recur. The same handful of problems account for most of what goes wrong
in a batch, they are invisible until something downstream fails in an unrelated-looking way, and
without a record each one gets rediagnosed from scratch the next time. A row here is also the
raw material for a report to the archive — see [Reporting upstream](reporting-upstream.md).

**Every row is a measurement, not an impression.** The `evidence` column names a file, a work
directory, or a published [post-mortem](post-mortems.md) that still contains the error text or
the numbers quoted. If you cannot point at
something, the row does not go in. This is also what makes a row usable in a helpdesk ticket:
archives act on reproducible specifics and close vague reports.

The two tables below are separated by a distinction worth keeping sharp:

- **[Defects](#defects)** — the archive's data or metadata is wrong, incomplete, or internally
  inconsistent. A depositor or the archive has to fix it. Worth reporting.
- **[Correctly labelled, rejected late](#correctly-labelled-rejected-late)** — the submission is
  accurate and says plainly what it is. *We* downloaded it anyway and rejected it afterwards.
  Nothing to report; this is a pipeline gap, tracked in [Failure modes](failure-modes.md).

Mixing the two is the most expensive mistake available here. It produces helpdesk tickets that
are wrong on their face and hides real pipeline bugs behind "bad public data".

## Defects

| Accession | Archive | Pathology | What was measured | Detected by | First seen | Evidence | Reported |
|---|---|---|---|---|---|---|---|
| GSE247111 | GEO / SRA | `shared-biosample` | Each GEX library and its CellPlex (CMO) tag library share one BioSample, so two GSMs resolve to an identical SRS. A run→sample mapping keyed on BioSample silently merges a tag library into a GEX sample. R2 complexity over the first 50,000 reads separates them: SRR26669510 has 651 distinct 30-mers, top one 82.8% (tag); SRR26669511 has 15,170, top one 5.9% (GEX). Merged outputs show 32–39% valid barcodes against 97.9% for a healthy sample in the same run. | R2 complexity probe; valid-barcode rate in the Cell Ranger summary | REQ-74217 / run 337224 | [Run 337224 post-mortem](https://claude.ai/code/artifact/7d89d7e0-705a-4c5b-bdb3-621615a9504c); archived run record for that run | no |
| GSE200629 | SRA | `zero-length-read` | `fastq-dump` discards an entire read as zero length in the archive, across **24 runs**. Verbatim: `ERROR: fastq-dump discarded a whole read of SRR18718349 for having zero length in the archive.` | `SRA2FASTQ` | batch 5 / `tender_brattain` | `failed5.log` names the runs, not the series: grep `SRR18718349`. Series via `searchlist.tsv` (SRR->GSM6040360) then `allhumandatasets.tsv` (GSM->GSE) | no |
| GSE202052 | SRA | `zero-length-read` | Same signature as GSE200629, **8 runs**; 32 across the two series in one batch. | `SRA2FASTQ` | batch 5 / `tender_brattain` | `failed5.log`, e.g. `SRR19039063` (GSM6090510); same two-step mapping as above | no |
| GSE206528 | GEO / SRA | `no-run-id` | `GSM6255907` is present in the series but has no experiment or run ID in the SRA export, so the sample cannot be resolved to anything downloadable. Verbatim: `ERROR: No experiment or run ID found for GSM6255907 in GSE206528.sra.tsv!` | `FETCH10XMETA` | batch 6 / `spontaneous_ampere` | `failures6.tsv`; `nf-work/c6/3b76ea1d9af7fb58568cdd1ea39d5a` | no |
| E-MTAB-10581 | ArrayExpress | `missing-files` | 72 FASTQ URIs named in the SDRF do not resolve. All 24 runs of the dataset failed at download; the dataset produced no output at all. | `WGET10X` | REQ-74217 | [REQ-74217 post-mortem](https://claude.ai/code/artifact/bcb35389-362e-4e89-9309-597e073aa66b) | no |
| E-MTAB-8060 | ArrayExpress / BioStudies | `listing-disagrees` | The BioStudies file listing and the SDRF disagree about which FASTQs exist, so a validation trusting BioStudies alone concludes the FASTQs are absent and reroutes the dataset. 15 runs affected. Part cause is ours (the URI check added in `3e5ec78`), but the two archive sources genuinely differ. | ArrayExpress URI validation | run 337224 | [Run 337224 post-mortem](https://claude.ai/code/artifact/7d89d7e0-705a-4c5b-bdb3-621615a9504c) | no |
| SRR19391245 | SRA | `single-mate-dump` | `fastq-dump` returns one FASTQ where the run is paired and two or more are expected: `ERROR: fastq-dump produced 1 FASTQ file(s) for SRR19391245`. The August 2026 run saw this at scale — 89 runs `single-fastq-from-sra`, plus 48 `read-counts-differ` where mates dumped at different depths. | `SRA2FASTQ` | batch 6 / `spontaneous_ampere`; at volume August 2026 | `failures6.tsv`; `failed.log` (August) | no |
| PRJCA017779 | GSA | `unsupported-archive` | A GSA (Genome Sequence Archive, China) accession with no ENA mirror, so no route resolves it. Not a defect in the submission — a gap in archive coverage — but it fails identically and is worth recording so it is not rediagnosed. | `FETCH10XMETA` | REQ-74217 | [REQ-74217 post-mortem](https://claude.ai/code/artifact/bcb35389-362e-4e89-9309-597e073aa66b) | n/a |

## Correctly labelled, rejected late

These submissions are accurate. `library_strategy` says what they are, and the pipeline
downloaded them before checking. **Do not report these.** The cost column is the argument for
the pre-download screen tracked in [Failure modes](failure-modes.md).

| Accession | Archive | What it actually is | Labelled in the archive as | Runs | Downloaded |
|---|---|---|---|---|---|
| GSE182791 | GEO | Bulk RNA-seq, 2×151 and 2×101. A SuperSeries: 5 of 75 samples are 10x, the other 70 are bulk replicates | `library_strategy = RNA-Seq`, no Cell Ranger in the processing description | 70 | 334.2 GB |
| GSE208195 | GEO | Bisulfite-seq (snmC), 2×150 | `library_strategy = Bisulfite-Seq` | 7 | 173.6 GB |
| GSE239932 | GEO | scATAC-seq half of a Multiome study, 50+49 | `library_strategy = ATAC-seq` | 8 | 60.7 GB |
| GSE218314 | GEO | scATAC-seq half of a Multiome study, 50+49 | `library_strategy = ATAC-seq` | 6 | 23.5 GB |
| GSE191286 | GEO | scATAC-seq, 2×50 | `library_strategy = ATAC-seq` | 4 | 4.8 GB |
| GSE109816 | GEO | Not 10x at all; no barcode read, best whitelist match 0.6% across all eight whitelists and both offsets | — | 880 | — |

596.8 GB across the first five, in a single run, to reject data whose own metadata said it was
not gene expression. GSE109816 alone was 880 of the 1814 failures in August 2026.

> **The GEX halves of GSE218314 and GSE239932 are fine** and should be processed. Only the ATAC
> halves belong in this table. Both were nonetheless lost to a separate, unrelated bug — Cell
> Ranger refusing auto-detected ARC-v1 — which is a pipeline problem, not an archive one.

## Adding a row

1. Establish the fact from the run's own artifacts — the archived `failed<N>.log` block, the
   `failures<N>.tsv` row, or a measurement taken in the work directory. Quote the error
   verbatim; paraphrase loses the string the next person will grep for.
2. Decide which table it belongs in by asking one question: **would the depositor or the archive
   have to change something?** If the submission is accurate and we simply did not read it, it
   is not a defect.
3. Fill `evidence` with something that still exists. A `nf-work` hash is fine while the work
   directory survives; once a run is archived, prefer its path under
   `/nfs/cellgeni/reprocessing-runs/`.
4. Leave `Reported` at `no` until a report is actually filed, then replace it with the link.

Past roughly 50 rows this should become a TSV with its own ignore exception; until then the
table stays readable in a diff and on GitHub.
