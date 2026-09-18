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
| GSE212964 / GSE212965 | GEO / SRA | `split-mate-runs` | R1 and R2 of one sequencing run are registered as **two separate SRA run accessions**, for 5 samples (`GSM6566486`-`GSM6566490`; run pairs SRR21492342/343, 371/372, 373/374, 375/376, 377/378). Each accession holds two 10 bp index reads plus one 101 bp read. In the first of each pair that 101 bp read matches the 3' v3 barcode whitelist at offset 0 at 94.3% (188542 of 199957 sampled); in the second it matches nothing at any offset. The halves are the same run: read names correspond record for record at the head of each file, same flowcell, lane, tile and cluster (`A00351:649:HFW57DSX3:1:1101:1434:1000`). Layout inference runs per *run*, upstream of the sample merge, so neither half can present a complete layout and both are rejected. The cDNA half additionally carries a uniform `?` quality string for every base. All 5 samples are recoverable; none was aligned. | `RENAME10XRUN` exit 2, as **complementary** errors — "Expected exactly one GEX biological R2 read, found 0" on the barcode half, "no whitelist hits" on the cDNA half | batch 8 / `cheeky_cuvier` | `failures8.tsv`; `nf-work/df/d31492a4d828435a231e8cb4d87a27` (SRR21492342) and `nf-work/20/a5d396b8444139697c3e08277d5ed6` (SRR21492343); `issues/2026-09-16-batch8-split-mate-runs.md` | no |
| GSE202476 | GEO / SRA | `split-mate-runs` | Same pathology, found retrospectively: `GSM6122993` has 8 failed runs in 4 strictly alternating complementary pairs (SRR19139475/476, 477/478, 479/480, 481/482). **Log text only** — batch 5's work dirs were deleted 2026-09-14, so the read-name and whitelist measurements that confirm the GSE212964 pairs cannot be repeated here. One sample lost. The other 3 rows of the same error in that batch (`SRR19096070/71/72`, GSE202308) have no complementary partner and are what a genuine ADT/HTO library looks like. | `RENAME10XRUN` exit 2, complementary-error signature applied to `failed5.log` | batch 5 / `tender_brattain` (recognised 2026-09-16) | `/nfs/cellgeni/reprocessing-runs/batch5/tender_brattain/` (archived `failed5.log`); `issues/2026-09-16-batch8-triage-verdict-split-mate.md` | no |
| GSE213370 / GSE213372 | SRA | `empty-run-object` | `GSM6584393` / `SRR21579492`: the SRA metadata declares a 2-read run of 100 + 50 bp, but the archived object is 63 KB and the dump yields a single 50 bp FASTQ of 17 reads (3.1 KB), with `dump.log` reporting "Read 1 spots" repeatedly. Verbatim: `ERROR: fastq-dump produced 1 FASTQ file(s) for SRR21579492: SRR21579492_1.fastq.` The guard fired correctly; the run holds no usable 10x data whatever the metadata says. Distinct from `single-mate-dump` in that the object is near-empty rather than merely short a mate. | `SRA2FASTQ` exit 1 | batch 8 / `cheeky_cuvier` | `failures8.tsv`; `nf-work/a0/b72ae4dc650b9960549e596bcd8962` (`dump.log`, the 17-read FASTQ) | no |
| GSE214914 | GEO / SRA | `unidentified-library` | Four samples (`GSM6617853`-`GSM6617856`) whose GEO titles carry no library label (`145790 AML`, `N34`, …) and whose reads are not gene expression. R1 151 bp matches the 3' v3 whitelist at 96.1%, so inference accepts the run; the 151 bp mate then aligns at **0.04% uniquely mapped, 99.92% "unmapped: too short"** against human 2020A, `exon_u` 0.00000, median 1-2 features per cell. What the library actually is was not identified. | `REPROCESS10X_STARSOLO10X` completing with an empty matrix; no task failed | batch 9 / `gigantic_hypatia` | `results/batch9/mapping_qc_stats.tsv`; `nf-work/cf/c03024faa221ee3fc7dff91d753de8` (`GSM6617853`, incl. `Log.final.out`) | no |
| GSE216883 | GEO | `mislabelled-library` | `!Sample_title` says `scRNAseq` for three samples (`SFMC, AS, scRNAseq [SF1_S1]`, `[SF2_S2]`, `PB_S3`) whose geometry is a feature-barcode library: R1 26 bp matching 3' v2 at 97.6%, R2 **25 bp**. The dataset's genuine GEX samples are titled `[S1_GEX]`-`[S4_GEX]` and processed normally. Correctly rejected, but only the suffix distinguishes them — the word `scRNAseq` in the title is wrong. | `RENAME10XRUN` exit 2, `Expected exactly one GEX biological R2 read, found 0: []` | batch 9 / `gigantic_hypatia` | `failures9.tsv`; `nf-work/e1/8f8e94a8bae4e555946aed13206969` (SRR22105506) | no |
| GSE215120 | GEO / SRA | `heterogeneous-run-set` | Eight GSMs each bundle a genuine 10x run with one or two runs that are not 10x at all — 150 bp + 150 bp, no barcode read, best whitelist match 0.2%. Six further runs (`SRR21849467`-`472`) fail SRA extraction outright: `ERROR: fastq-dump discarded a whole read of SRR21849467 for having zero length in the archive.` The 10x run alone is aligned and the output is healthy (`GSM6622292`: 95.6% mapped, 2251 median features), so the sample survives — but between a third and two thirds of each sample's registered runs are silently dropped, with nothing in the output saying so. | `RENAME10XRUN` exit 2 (no whitelist hits) and `SRA2FASTQ` exit 1 | batch 9 / `gigantic_hypatia` | `failures9.tsv`; `nf-work/1a/73ee19a91b14fe400e495609f96259` (SRR21849525, rejected) vs the accepted `SRR20791206` at 97.3% v2 | no |
| GSE215915 | SRA | `near-empty-run` | Five runs whose FASTQ objects hold almost no data: sampling exhausts the file at 54, 65, 110, 264 and 19 records respectively, against a 200,000-record target. Geometry is 10 bp + 90 bp — no barcode read present. The affected samples still align from their three remaining runs. | `RENAME10XRUN` exit 2, `No supported 10x run layout matched whitelist evidence plus read geometry (no whitelist hits)` | batch 9 / `gigantic_hypatia` | `failures9.tsv`; `nf-work/fe/16613adfcd58227fc622e2b793a7c6` (SRR21933465) | no |
| GSE216673 | SRA | `single-mate-dump` | Five runs (`SRR22063067`-`071`), the whole dataset, dump one mate where two are expected: `ERROR: fastq-dump produced 1 FASTQ file(s) for SRR22063067: SRR22063067_1.fastq.` Same signature as SRR19391245 above. All 5 samples lost; the dataset produced no output. | `SRA2FASTQ` exit 1 | batch 9 / `gigantic_hypatia` | `failures9.tsv`, `failed9.log` | no |
| GSE215908 | GEO | `suspect-species` (**unconfirmed**) | All 16 samples are declared `Homo sapiens` (`!Sample_taxid_ch1 = 9606`) and align to human 2020A at only 21-29%, `exon_u` 0.11-0.14 — yet with 2700-11500 cells and 500-2300 median features, so the data itself is real and the matrices are usable, not empty. Sample titles describe transplanted and co-transplanted tumour models. A xenograft or a mislabelled organism would both look like this; **neither was confirmed**, and no mouse alignment was run to test it. Recorded so the next batch does not rediagnose it from scratch. | `results/batch9/mapping_qc_stats.tsv`; no task failed | batch 9 / `gigantic_hypatia` | `mapping_qc_stats.tsv` rows `GSM6645415`-`GSM6645430` | no |

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

### Accurately labelled, and not rejected at all

A third case, found in batch 9 and worse than rejecting late: these submissions say what they
are, the pipeline accepted them as gene expression anyway, aligned them, and **published the
resulting near-empty matrices** alongside real data. Nothing failed, so none of it appears in a
failure list.

The library type is stated in `!Sample_title` and nowhere machine-readable: all 38 samples carry
`!Sample_library_strategy = RNA-Seq` and `!Sample_library_source = transcriptomic single cell`,
which for a hashing or antibody library is not even untrue. Read geometry does not separate them
either — in GSE216914 the CSP and VDJ libraries of one sample share an identical
`10 / 10 / 26 / 90` layout, and one maps at 0.9% while the other maps at 93.9%.

| Accession | What it actually is | Titled in GEO as | Samples | `all_u+m` | `exon_u` |
|---|---|---|---|---|---|
| GSE216999 | Cell hashing (HTO) | `Library N cell hashing` | 19 | 0.008-0.388 | 0.0004-0.0026 |
| GSE216914 | Cell-surface protein (ADT/CSP) | `PSA1 CSP`, `PSA2 CSP`, … | 4 | 0.005-0.009 | 0.00000 |
| GSE216858 / GSE216859 | Feature barcode, unspecified | `Custom library_48h`, `_1wk`, `_2wk` | 3 | 0.019-0.044 | 0.014-0.040 |
| GSE215842 / GSE217326 | Hashtags | `scrnaseq with hashtags … [HT]` | 2 | 0.007-0.034 | 0.002-0.021 |
| GSE215253 | CRISPR guide enrichment | `… enrichment PCR` | 2 | 0.006-0.010 | 0.0000045 |
| GSE217499 | Hashtags | `… hashtag library` | 2 | 0.011-0.012 | 0.0006-0.0007 |
| GSE216040 | CRISPR guide enrichment | `gRNA enrichment library` | 1 | 0.009 | 0.0007 |
| GSE215825 | Not established | `d40, scRNAseq,pooled` | 1 | 0.002 | 0.00000 |

Real samples in the same batch run 0.66-0.98 `all_u+m` and 0.43-0.72 `exon_u`. Adding the four
`unidentified-library` samples of GSE214914 from the Defects table gives 38 of batch 9's 784
aligned samples — 4.8% — published as matrices of essentially nothing.

The rest of GSE215253 is the counter-example worth keeping in view: its 54 enrichment-PCR runs
*were* correctly rejected at `RENAME10XRUN`, and its 41 scRNA-seq samples all processed cleanly.
Only two of its feature libraries slipped through, because their biological read is 42 bp rather
than the 20 bp of the rest. No length threshold separates that from good data —
`GSM6660161` (56 bp), `GSM6668926`, `GSM6681083` and `GSM6681046` (55 bp) all map at 78-94%.

Tracked, with the measured gate that would catch them, in
`issues/2026-09-18-batch9-no-mapping-rate-gate.md`. Nothing here is worth reporting to GEO.

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
