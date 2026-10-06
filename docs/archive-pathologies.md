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
| GSE247111 | GEO / SRA | `shared-biosample` | Each GEX library and its CellPlex (CMO) tag library share one BioSample, so two GSMs resolve to an identical SRS. A run→sample mapping keyed on BioSample silently merges a tag library into a GEX sample. R2 complexity over the first 50,000 reads separates them: SRR26669510 has 651 distinct 30-mers, top one 82.8% (tag); SRR26669511 has 15,170, top one 5.9% (GEX). Merged outputs show 32–39% valid barcodes against 97.9% for a healthy sample in the same run.  **Recurred in batch 21 (`sick_galileo`) through `FETCH10XMETA`'s BioSample fallback**: 12 of 22 requested GSMs lost at exit 0, and all 10 published matrices carry a `CellPlex files` GSM's name while pooling that condition's GEX run. | R2 complexity probe; valid-barcode rate in the Cell Ranger summary; requested-vs-emitted check | REQ-74217 / run 337224 | [Run 337224 post-mortem](https://claude.ai/code/artifact/7d89d7e0-705a-4c5b-bdb3-621615a9504c); archived run record for that run | no |
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
| GSE218936 / GSE218937 | GEO / SRA | `shared-run-set` | All 8 GSMs of the series resolve to just **two** sets of four SRA runs — `GSM6760291`-`294` to `SRR22356862/863/874/875`, `GSM6760295`-`298` to `SRR22356858/859/860/861` — although the eight are distinct experimental conditions (1G vs µG, stimulated vs unstimulated, two donors). The pipeline emitted only the last GSM of each group and aligned it against all four runs, so two published matrices each pool four conditions under one sample's name, at 96% mapping and 1865 median features. **Resolved 2026-09-27:** GEO gives each of the eight GSMs its own SRX, and each group of four shares one BioSample (`SAMN31816603`, `SAMN31816604`). The pooling is ours: `collect_metadata.sh` falls back to a BioSample lookup when the GSM is absent from the SRA table. The conditions are separate libraries and can be recovered. See GSE234714 below. Same family as GSE247111 above; different consequence. | requested-vs-emitted sample check (`SKILL.md §4`); `sample_x_run.tsv` in the run's own metadata | batch 10 / `trusting_carson` | `results/batch10/metadata/GSE218936/{sample_x_run.tsv,links.tsv}`; `nf-work/14/f604b89d3f82ab48ee40ffd3f622bd`; `issues/2026-09-19-batch10-metadata-drops-samples-silently.md` | no |
| GSE218316 | GEO / SRA | `unidentified-library` | All 5 samples. R1 151 bp matches the 3' v2/v3 whitelist at 95.6–97.2% so inference accepts the run, and the 151 bp mate then aligns at **0.76% uniquely mapped, 99.09% "unmapped: too short"** against human 2020A, `exon_u` 0.0002–0.003, 3–12 median features. Declared `Homo sapiens`, titles are plain `Donor N - Nh treatment` with no library-type marker, `library_strategy = RNA-Seq`. The dataset produced no usable output at all. What the library actually is was **not identified** — same signature as GSE214914 above. | `REPROCESS10X_STARSOLO10X` completing with an empty matrix; no task failed | batch 10 / `trusting_carson` | `results/batch10/mapping_qc_stats.tsv`; `nf-work/82/f17836b943708b30ab8610ec5c1e2e` (`GSM6740283`, incl. `Log.final.out`) | no |

| GSE220774 | GEO / SRA | `concatenated-read` | Six samples (`GSM6817217`/`218`/`219`/`222`/`223`/`224`, one run each: `SRR22696361`/`363`/`364`/`366`/`368`/`369`) were deposited with R1 and R2 **joined into a single 150 bp record**, and the 8 bp sample index uploaded as the `_2` file. The first record of `SRR22696361_1.fastq.gz` reads `ANAGAGCCACCGTACG` (16 bp CB) + `CCTCCCTGCGAG` (12 bp UMI) + a 30 bp polyT + cDNA; `_2` is `AGATCTCG` on every record. Whitelist match on the first 16 bp is 94.2% against 3' v3, so the barcode is real and the data is recoverable in principle by splitting at base 28 — but no biological read exists as deposited. ENA lists exactly two files per run, so nothing was missed at download. The dataset's other 9 samples processed normally. 6 samples lost. | `RENAME10XRUN` exit 2, `Expected exactly one GEX biological R2 read, found 0: []` — **note this is the same message three unrelated causes emit**, see `issues/2026-09-16-batch8-triage-verdict-split-mate.md` | batch 11 / `irreverent_elion` | `failures11.tsv`; `nf-work/d4/c1248d16b75341af8a983dc2c4ca3c` (SRR22696361, incl. `.command.err` and the staged FASTQs) | no |

| GSE229166 | GEO / SRA | `no-run-id` | `GSM7156182` has no experiment or run ID in either the SRA or the ENA export. Same archive defect as GSE206528 above, but recorded separately because the **consequence is the whole dataset**: `collect_metadata.sh` `break`s out of its per-sample loop on the first unresolvable GSM, so the other 20 samples were never looked up and the series produced no `links.tsv` at all. 21 samples lost — 75% of batch 14's entire sample loss, from one bad record. Verbatim: `ERROR: No experiment or run ID found for GSM7156182 in GSE229166.sra.tsv!` then `…ena.tsv!` then `ERROR: Failed to get sample, experiment, and run IDs for GSE229166 using any of the available metadata files!` | `FETCH10XMETA` exit 1, 5 attempts | batch 14 / `maniac_keller` | `failures14.tsv`; `nf-work/c7/a87a63f6ab9f91bc17d84a9268e3eb`; `issues/2026-09-25-batch14-fetch10xmeta-aborts-whole-dataset.md` | no |
| GSE228428 | ENA | `ena-fastq-omits-barcode-read` (**re-diagnosed 2026-09-30**, first filed as `barcode-read-absent`) | Five runs (`SRR23998252`/`254`/`256`/`258`/`260`), 5 samples, the whole dataset. ENA's two-file export holds a 10 bp `_1` (the sample index, constant across every record: `AGACCATCGG` in all of `SRR23998252`) and a 51 bp `_2`, so there is nothing for a whitelist to hit. **The barcode read was submitted and SRA has it.** SRA's run record lists the submitter's `CR75-5-scR_S5_L002_{I1,I2,R1,R2}_001.fastq.gz` and stores four reads per spot at 10, 10, **28**, 51 bp. ENA's export dropped the 28 bp barcode read. First filed here as a submitter defect, "not recoverable as deposited". That was wrong: the 5 samples are recoverable from SRA. Same defect as GSE241998/GSE241999 below. | `RENAME10XRUN` exit 2, `No supported 10x run layout matched whitelist evidence plus read geometry (no whitelist hits)`; SRA `sra-db-be/run_new?acc=SRR23998252` for the read structure | batch 14 / `maniac_keller` | `failures14.tsv`; `nf-work/81/f9f7d575ecc8ead7050eb49ffa60ce`; `issues/2026-09-30-batch19-ena-fastq-omits-barcode-read.md` | no |
| GSE232523 | SRA | `near-empty-run` | Five runs (`SRR24560599` and four siblings), the whole dataset. Each run's three FASTQ objects total 27-33 KB and hold **318 records**, against a 200,000-record sampling target. `gzip -t` passes on all of them, so this is the deposit, not a truncated download. The geometry is right — `_1` 8 bp index, `_2` 28 bp barcode, `_3` 91 bp cDNA — but 316 usable barcodes give a best whitelist hit of 16.1% (v2), far below threshold and statistically meaningless at that depth. A correct rejection that `triage.py` nonetheless files under `chem-no-whitelist-hit` / `investigate`, identically to GSE228428 above, because the record count is not in `first_error`. | `RENAME10XRUN` exit 2, `No supported 10x run layout matched whitelist evidence plus read geometry (no whitelist hits)` | batch 15 / `infallible_dalembert` | `failures15.tsv`; `nf-work/d8/c495f36f5a613b881f0307f8fd4e75` | no |
| GSE227640 | SRA | `unmatched-barcode-set` (**unconfirmed**) | Two runs (`SRR23904484`, `SRR23904496`), one of the two runs of each of `GSM7104497` and `GSM7104492`. Geometry looks like a clean 10x v2 library — `_1` 26 bp (16 CB + 10 UMI), `_2` 98 bp — and 116,777 barcodes were usable, so depth is not the problem. But the best whitelist hit is **19.6%** against `gex_3pv2_or_5pv1v2`, with every other chemistry at 0.0-2.2%. That is far above noise and far below the ~90% a real match gives, and **what the barcode set actually is was not determined** — a custom or degenerate barcode set, or a non-10x droplet platform, would both look like this. **No samples were lost**: both GSMs aligned from their surviving run (`exon_u` 0.55 and 0.67, 2,478 and 1,016 median features), which is also why this is easy to miss. Recorded so it is not rediagnosed from scratch. | `RENAME10XRUN` exit 2, `No supported 10x run layout matched whitelist evidence plus read geometry (no whitelist hits)` | batch 14 / `maniac_keller` | `failures14.tsv`; `nf-work/fd/f7405f80c0d79a9239e6f574ce6832` | no |
| GSE229170 | SRA | `zero-length-read` | Four runs (`SRR24108136`-`139`), all of `GSM7156264`, which was the only sample lost. `fastq-dump` discards an entire read as zero length. Verbatim: `ERROR: fastq-dump discarded a whole read of SRR24108136 for having zero length in the archive.` Same signature as GSE200629 and GSE202052 above. | `SRA2FASTQ` exit 1 | batch 14 / `maniac_keller` | `failures14.tsv`, `failed14.log` | no |
| GSE229279 | GEO / SRA | `library-type-mismatch` (**unconfirmed**) | 12 samples, paired as 6 `…, transcriptome` and 6 `…, TCR`. Neither half behaves as titled. The six `TCR` samples carry **no V(D)J enrichment at all** — TR/IG segment UMI fraction 0.04-0.67% against a measured GEX baseline of 0.28% — yet map poorly for GEX (`exon_u` 0.065-0.26, 329-863 median features). Meanwhile `GSM7157924`, titled `Healthy control skin patient 2, transcriptome`, aligned to nothing (`exon_u` 5.8e-07, 1 median feature, 1,465 cells) although its five `transcriptome` siblings are healthy (`exon_u` 0.33-0.60, 665-1,808 features). A swapped GSM→run relation would explain part of it but not the absent enrichment; **no cause was established**. No task failed, so none of this appears in any failure list. | matrix V(D)J-segment fraction; `results/batch14/mapping_qc_stats.tsv` | batch 14 / `maniac_keller` | `results/batch14/starsolo/GSE229279/*/output/Gene/filtered`; `results/batch14/mapping_qc_stats.tsv` | no |
| GSE234714 / GSE234717 | GEO / SRA | `shared-biosample` | 20 GSMs (`GSM7473691`-`GSM7473710`, four replicates for each of five conditions) share **five BioSamples, four GSMs each**. Each GSM has its own SRX and each SRX has exactly one run, so the archive is unambiguous. The SRA sample name is the submitter's (`1023ZDH`, `D2F4`, …), never the GSM. `FETCH10XMETA` therefore resolved every GSM by BioSample, gave each one all four of its group's runs, and emitted only the last GSM of each group. 15 samples lost with exit 0, and **5 published matrices each pool four replicate libraries**. The SRX lookup that fixes it is in `issues/2026-09-19-batch10-metadata-drops-samples-silently.diff`. GSE234717 is the duplicate row. Its own BioProject holds none of these runs, and it failed loudly with `No experiment or run ID found for GSM7473691`. A shared BioSample is legal in SRA, so this is a defect in our pipeline, recorded here because of how the submission is shaped. | requested-vs-emitted sample check; `sample.relation.list` against `sra.tsv` | batch 16 / `sad_brenner` | `results/batch16/metadata/GSE234714/{sample.relation.list,sra.tsv,sample_x_run.tsv,links.tsv}`; `nf-work/0d/040d1d558387224480e40a7a4ac8b1` (GSE234717) | n/a |
| GSE233338 / GSE233339 | SRA | `index-reads-only` | 4 samples × 3 runs (`GSM7425875`-`878`, the GEX arm of a Multiome study), the whole dataset. Every run holds **two 10 bp reads and nothing else**. ENA's file report for `SRR24727163` gives 31,095,121 spots and 621,902,420 bases, i.e. 20 bp per spot. Only the i7/i5 sample indices were submitted; there is no barcode read and no cDNA. 0 whitelist hits at any offset. Not recoverable as deposited. Same family as GSE228428 above, but there the cDNA was present. | `RENAME10XRUN` exit 2, `no whitelist hits`, `len=10` on both files | batch 16 / `sad_brenner` | `failures16.tsv`; `nf-work/0e/028c9bc25092d8510d03e67d79a292` (SRR24727163) | no |
| GSE235534 | SRA | `index-reads-only` | 12 samples, one run each, the whole dataset (3' v3.1 GEX per the protocol text). Same shape as GSE233338: `SRR24988626` is 108,762,547 spots × 20 bp, `_1` and `_2` both 10 bp. 12 samples lost. | `RENAME10XRUN` exit 2, `no whitelist hits` | batch 17 / `cranky_cuvier` | `failures17.tsv`; `nf-work/d1/f9467c4cf3a2b5a48c15b805ed3ddd` (SRR24988626) | no |
| GSE235330 / GSE235332 | SRA | `index-reads-only` | One of the two runs of each of `GSM7500564`, `GSM7500569` and `GSM7500573` (`SRR24975949`, `…940`, `…932`) is registered `SINGLE` with **8 bases per spot**: 289,661,196 spots and 2,317,289,568 bases for `SRR24975949`. That is the i7 index alone. `fastq-dump` correctly yields one FASTQ. All three samples were published from their other run, so **half of each sample's sequencing is missing** from its matrix and nothing records that. | `SRA2FASTQ` exit 1, `fastq-dump produced 1 FASTQ file(s)` | batch 17 / `cranky_cuvier` | `failures17.tsv`; `nf-work/3a/14a77018d82706597a3b7533874a65` (SRR24975949) | no |
| GSE236107 | GEO / SRA | `split-mate-runs` | `GSM7518380`'s runs `SRR25069197` and `SRR25069198` each hold **75,704,350 spots**. `…197` carries the 28 bp barcode read (3' v3 at 96.0%) plus an 8 bp index. `…198` carries an 8 bp index plus the 101 bp cDNA. Same pathology as GSE212964 above, with the complementary-error signature. This time the sample's other 6 runs were complete, so it was published from 6 of 8 runs rather than lost. | `RENAME10XRUN` exit 2, complementary errors | batch 17 / `cranky_cuvier` | `failures17.tsv`; `nf-work/87/0782e0e7edc9dcfb7057a5085ea4fe` (SRR25069197); `issues/2026-09-16-batch8-split-mate-runs.md` | no |
| GSE234963 / GSE234971 | SRA | `partial-mate` | `SRR24937677` (`GSM7486861`, the sample's only run) declares 1,051,028,696 spots, but its third read dumps **128,188,908 records**, 12% of the spot count. Consistent across 10 attempts, including the duplicate row, so it is the object and not a transient. Verbatim: `ERROR: mates of SRR24937677 disagree on record count (SRR24937677_3.fastq has 128188908, expected 1051028696).` 1 sample lost. | `SRA2FASTQ` exit 1 | batch 16 / `sad_brenner` | `failures16.tsv`; `nf-work/e6/43cfd05df062ae8e8c71aa5a62f24b` | no |
| GSE241998 / GSE241999 | ENA | `ena-fastq-omits-barcode-read` | 3 of 13 samples (`GSM7747454`, `GSM7747457`, `GSM7747458`; runs `SRR25820081`, `SRR25819676`, `SRR25819675`). SRA stores each spot as 3 or 4 reads, cDNA 91 bp, barcode 28 bp and one or two index reads, **in an order that varies from run to run**, and holds the submitter's original `I1`/`R1`/`R2` files. ENA exports two files. Where the barcode read sat between the others in SRA's order (`SRR25819675`: 91, **28**, 8), ENA served 91 + 8 and dropped it. Where it sat last (`SRR25819677`: 91, 8, **28**), ENA served 91 + 28 and the run aligned. The 10 siblings came out fine by luck of read order. The pipeline prefers ENA whenever a `_1`/`_2` pair exists, so the 3 samples were lost although SRA has everything. Report to ENA: its export differs from the SRA object. | `RENAME10XRUN` exit 2, `no whitelist hits`, with files at `len=91` and `len=8`/`10`; SRA `sra-db-be/run_new` for read order | batch 19 / `focused_goldberg` | `failures19.tsv`; `nf-work/83/6b6726be79296bc453b0810a98ce3c` (SRR25819675), `nf-work/78/caa683aeb9f32a7de8ef8c57a6ee14` (SRR25820081); `issues/2026-09-30-batch19-ena-fastq-omits-barcode-read.md` | no |
| GSE241739 / GSE242039 | SRA | `index-reads-only` | Each of the 6 GSMs (`GSM7734667`-`672`) has 4 runs: two real libraries at 118 bp per spot (28 + 90) and, **registered as separate runs, their index reads** at 20 bp per spot (10 + 10) with the identical spot count (`SRR25775008` and `SRR25775009` both hold 176,038,835 spots). 12 index-only runs rejected per batch, 24 across the two. **Harmless**: every sample was aligned from its two real runs, and the partial-input check flags all 6 only because the runs it counts include the index-only ones. GSE242039 is the same series under a second accession, processed in batch 19. | `RENAME10XRUN` exit 2, `no whitelist hits`, both files `len=10` | batch 18 / `compassionate_shaw` | `failures18.tsv`, `failures19.tsv`; `nf-work/cc/716c0ca779a8d1a1f87990674a0450` (SRR25775007); `results/batch18/metadata/GSE241739/GSE241739.sra.tsv` | no |
| GSE241292 | GEO / SRA | `shared-biosample` + `unloaded-run` | 13 requested GSMs, each a "main sample" (`MS`, short-read) or a "sub-sample" (`SS`, long-read) library. Each MS/SS pair shares **one BioSample**, and each GSM has its own SRX. `FETCH10XMETA` resolved by BioSample, gave both GSMs of a pair both runs, and emitted only one. **5 MS samples lost at exit 0** (`GSM7720737`/`738`/`739`/`745`/`746`), and **5 published SS matrices each pool the MS library's reads with their own** (`GSM7720741`/`742`/`743`/`749`/`750`; `GSM7720741` = 161 M SS + 691 M MS reads). Separately, 3 MS runs (`SRR24186411`/`424`/`425`) are registered with **0 spots and 0 bytes** and no download URL in SRA or ENA. They failed as `wget-empty-url` on 5 attempts each. Their SS partners (`GSM7720744`/`751`/`752`) were published from their own run alone, which is correct by accident. | requested-vs-emitted sample check; `sample.relation.list` against `accessions.tsv`; `WGET10X` exit 1 `no download URL for SRR… (type SRA)` | batch 18 / `compassionate_shaw` | `results/batch18/metadata/GSE241292/{sample.relation.list,accessions.tsv,sra.tsv,links.tsv}`; `nf-work/3c/cb9c6d693c5d68a4f75f3641630fce` (SRR24186411) | n/a (shared BioSample: our defect); no (unloaded runs) |
| GSE233279 | GEO / SRA | `failed-library` (**unconfirmed**) | `GSM7040229` (`C-13291`, caudate snRNA-seq), one of 24 samples. Its matrix calls **676,629 cells at 2 median features**. The other 23 call 1,111-12,088 cells at 2,482-3,543 features, except `GSM7040228` at 81,811 / 430. The barcodes are valid (97.7% whitelist, 92.4% unique mapping to genome), but the reads spread thinly: 13 reads and 2 UMIs per barcode, saturation 0.81. The strand test chose `Reverse` (fwd 41%, rev 54%) where every sibling is `Forward`. A library with almost no cell-associated signal would look like this. Whether it is a failed library or something we did was **not established**. `exon_u` is 0.056, just above the 0.05 screen, so nothing flagged it. | `mapping_qc_stats.tsv` (Cells, Med_nFeature); `strand.txt` | batch 16 / `sad_brenner` | `results/batch16/starsolo/GSE233279/GSM7040229/{strand.txt,output/GeneFull/Summary.csv}`; `nf-work/05/6f3b8ab8a298a794f915434d5130dc` (renamed FASTQs) | no |
| GSE245998 | SRA | `split-read-runs` | Multiome snRNA-seq, GEX arm, 6 GSMs (`GSM7853050`-`055`). Each GSM's one SRX holds **four runs of one read each**, with identical spot counts and identical read names. For `SRX22172188` (`GSM7853055`) the runs are `SRR26468217` 10 bp, `SRR26468218` 10 bp, `SRR26468183` **30 bp**, `SRR26468184` **98 bp**, all 264,559,316 spots, all starting `A00351:600:HM3CMDSX2:3:1101:1850:1000`. Run aliases are `GSM7853055_r1`..`_r4`. So I1, I2, R1 and R2 were registered as separate runs. The 30 bp read matches `gex_737K-arc-v1.txt` at 94.1% (record 0) and 94.4% (record 2,000,000). ENA serves one unsuffixed file per run, so the pipeline takes the SRA route and every run fails `SRA2FASTQ`'s two-mate guard, which is correct for any single run. The four-way form of `split-mate-runs` below. Whole dataset lost (6 samples), recoverable by pairing r3 + r4 per SRX. 77.9 GB of SRA objects pulled, 5 attempts each. | `SRA2FASTQ` exit 1, `fastq-dump produced 1 FASTQ file(s)`; runs of one SRX with equal spot counts in `sra.tsv` | batch 20 / `sleepy_curie` | `failures20.tsv`; `nf-work/e4/5119f59f1169d8facf726999881af3` (SRR26468183, incl. `dump.log`); `results/batch20/metadata/GSE245998/GSE245998.sra.tsv`; `issues/2026-09-16-batch8-split-mate-runs.md` | no |
| GSE247205 | GEO / SRA | `split-mate-runs` | 9 GSMs (`GSM7885251`-`259`), 2 runs each, in one SRX with identical spot counts (`SRR26701901`/`902`: `SRX22401521`, 229,202,779 spots each). The first run of each pair holds a 10 bp index + a 101-107 bp barcode read matching 3' v3 at 95.6-95.7%; the second holds the same index + the cDNA read. Read names correspond (`A00674:613:HCHWJDSX7:4:1101:1063:1016` heads both). Same shape as GSE212964 above. Whole dataset lost, recoverable. | `RENAME10XRUN` exit 2, complementary errors: `Expected exactly one GEX biological R2 read, found 0` / `no whitelist hits`, 9 + 9 | batch 21 / `sick_galileo` | `failures21.tsv`; `nf-work/8f/4fc9dd97b5634d2fe7ef73a4eaf410` (SRR26701901), `nf-work/ea/b1a9e17e4e071b9d9dc2fe51708c9e` (SRR26701902) | no |
| GSE246613 | GEO / SRA | `shared-biosample` + `unloaded-run` | 251 requested GSMs in 104 BioSamples (57 × 2 GSMs, 45 × 3, 2 × 1): per patient and timepoint a `_P` GEX library, a `_TCR` V(D)J library and sometimes an `_N` GEX library, each with its own GSM and SRX, sharing one BioSample. `FETCH10XMETA` resolved by BioSample and emitted the last GSM of each group, which is the `_TCR` one in 101 of 102 multi-GSM groups. **147 samples lost at exit 0**; **102 published matrices pool 2-3 libraries** under the survivor's name (`GSM7872699` `h01A_TCR` = `SRR26541224` TCR, 7.6 M spots + `SRR26541225` `h01A_P` GEX, 87.4 M spots; 29.8 M UMIs, 3.7% in TR/IG segments). An SRX lookup resolves 251/251 to their own single run. Separately, `SRR26541168` (`h30C_P`, `GSM7872821`, not requested) has **0 spots** and no ENA file; its BioSample partner `GSM7872822` (`h30C_TCR`) was given it, failed it at `WGET10X`, and was published from its TCR run alone. | requested-vs-emitted check; `sample.relation.list` against `sra.tsv`; `WGET10X` exit 1 `no download URL for SRR26541168 (type SRA)` | batch 21 / `sick_galileo` | `results/batch21/metadata/GSE246613/{sample.relation.list,sra.tsv,sample_x_run.tsv,links.tsv}`; `issues/2026-09-19-batch10-metadata-drops-samples-silently.md` (update 2026-10-02) | n/a (shared BioSample: our defect); no (unloaded run) |
| GSE245310 | GEO | `lane-split-gsm` | 68 GSMs that are **14 libraries**: the submitter made one GSM per sequencing lane per index oligo (`DRG_GW8_sample1_lane1_S2-1` … `lane3_S2-4` is 12 GSMs of one library; `DRG8_GW14_S1_L001_R1_001` … `_004` is 4). Every published matrix is a fraction of a library, and sibling GSMs share cell barcodes. Not what cost samples here: 52 of the 68 were lost to the STARsolo strand-test segfault, which is our defect (`issues/2026-09-18-batch9-starsolo-strand-test-segfault.md`), and it split along library lines (7 libraries lost, 7 kept). | `!Sample_title` | batch 20 / `sleepy_curie` | `results/batch20/metadata/GSE245310/GSE245310_family.soft`; `failures20.tsv` | no |
| GSE249894 | ENA | `ena-fastq-corrupt-at-source` | Both FASTQs of `SRR27178495` (`GSM7966594`) and of `SRR27178635` (`GSM7966587`) are truncated gzip in ENA's store: `pigz -t` reports `incomplete deflate data` on every clean download. Each file matches ENA's own `fastq_bytes` **and `fastq_md5`** (e.g. `SRR27178495_1.fastq.gz`, 837,718,315 bytes, `7a4191c8cdfed9386ad3f34990689cb1`), so ENA checksummed the truncated objects. SRA lists both runs (34,829,780 and 9,129,703 spots × 310 bp); the SRA copies were not downloaded. Each sample was published from 11 of 12 runs (6.7% and 1.8% of reads missing). | `WGET10X` exit 1, 5 attempts each | batch 22 / `gloomy_church` | `failures22.tsv`; `nf-work/8e/a1be1e46afc976638733af4988652a`, `nf-work/6e/b317c9f15bbac0decee2c0955de617`, `nf-work/8e/26b14f7f75bb8201c976aff01009ff`, `nf-work/1c/a714458ca1ce2ac2312b1be1fd70e4`; ENA filereport queried 2026-10-04 | no |
| GSE249159 | SRA | `mixed-spot-layout` | `SRR27016858` (`GSM7927630`, 364,895,176 spots) and `SRR27016862` (`GSM7927627`, 295,018,806), each the sample's only run. Lane-1 spots hold 3 reads (8 bp index, 28 bp barcode, 91 bp cDNA), lane-2 spots 2 (28 bp barcode, 91 bp cDNA), so `--split-files` gives a `_1` mixing 8 and 28 bp, a `_2` mixing 28 and 91 bp, and a `_3` with half the spot count. Mean spot length 123 bp. Verbatim: `ERROR: mates of SRR27016858 disagree on record count (SRR27016858_3.fastq has 182634053, expected 364895176).` The data is complete; 2 samples lost. | `SRA2FASTQ` exit 1, 5 attempts each | batch 22 / `gloomy_church` | `failures22.tsv`; `nf-work/f3/f3eeedb3747a86e14725ab8d8e55de` (read heads per file) | no |
| GSE248489 | SRA | `mates-as-single-spots` + `mixed-spot-layout` | `GSM7915896`'s two runs, one SRX (`SRX22617800`). `SRR26923939` holds 405,579,164 **single-read** spots: the first half 28 bp, the rest 90 bp, with the same read names, i.e. R1 and R2 uploaded as two single-end files into one run (`fastq-dump produced 1 FASTQ file(s)`). `SRR26923940` (398,540,976 spots, mean 64 bp) dumps a `_1` alternating 10 bp and 90 bp records beside a 193,518,513-record 28 bp `_2`. Data complete; 1 sample lost. | `SRA2FASTQ` exit 1, 5 attempts each | batch 22 / `gloomy_church` | `failures22.tsv`; `nf-work/e3/b5a19f7bdb47394f371184f1564ad3`, `nf-work/83/d598b67aaae6f75375e76bdcb48860` | no |
| GSE254170 | GEO / SRA | `shared-biosample` | `GSM8427071` (`WaGa-RB1kd uninduced control`) and `GSM8427072` (`… dox-induced TA k.d.`) share BioSample `SAMN42510873`, with their own SRX each (`SRX25338934` → `SRR29841142`; `SRX25338935` → `SRR29841141`). `FETCH10XMETA` emitted only `GSM8427072` and gave it both runs: 1 sample lost at exit 0, and the published `GSM8427072` matrix pools control with knockdown. Our defect (GitHub #6), recorded here for the shared BioSample. | requested-vs-emitted check | batch 22 / `gloomy_church` | `results/batch22/metadata/GSE254170/` (`sample.relation.list`, `sra.tsv`, `links.tsv`) | no |
| GSE248788 / GSE249894 / GSE250444 | SRA | `index-reads-only` | Index reads registered as their own runs beside the library run, with the identical spot count in one SRX. GSE248788: one 20 bp-per-spot run (10 + 10) per GSM, 9 GSMs (`SRR26978583`-`599`, odd numbers), beside a 116 bp-per-spot run. GSE249894: 7 runs `SRR33413768`-`780` (even numbers) of `GSM8964093`-`096`, beside 300 bp-per-spot runs. GSE250444: `SRR27476635`, `SRR27476638`, 8 bp per spot. Harmless: every sample aligned from its real runs. | `RENAME10XRUN` exit 2 (`no whitelist hits`), `SRA2FASTQ` exit 1 (1 file) | batch 22 / `gloomy_church` | `failures22.tsv`; `GSE248788.sra.tsv`, `GSE249894.sra.tsv` | no |
| GSE248214 | SRA | `near-empty-run` | `SRR26881574`, one of `GSM7908456`'s 8 runs, holds 2,686 spots in SRA itself (8 + 28 + 100 bp, 0.2% best whitelist match). Harmless: 0.001% of the sample's 294 M spots. | `RENAME10XRUN` exit 2 | batch 22 / `gloomy_church` | `nf-work/58/284a8cbecf1c451230b6d8fa69fda5` | no |

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
| GSE233970 / GSE234040 | GEO | Smart-seq2, 2×148. One sample (`GSM7439998`, `SSC GrowthPlateSSC`), one run per cell | `!Sample_description = SmartSeq2`, protocol text `SmartSeq2 pipeline (Picelli et al.)` | 191 | not measured (1,528 tasks, ~10 CPU-h) |
| GSE245175 | GEO | SORT-seq (robotised CEL-seq2), plate-based, 28+90 on NextSeq 500 / NovaSeq. Not 10x: best whitelist hit 3.3% (3' v3) on the 28 bp read | `!Sample_extract_protocol_ch1 = We have made use of the SORT-seq protocol, which is a partially robotized version of the CEL-seq2 protocol`; titles `…, plate KK001` | 18 (9 samples requested) | 14.5 GB |
| GSE235787 | GEO | 10x Flex (Fixed RNA Profiling), 28+90. One sample, `GSM8086562`. Flex is rejected from the GEX path by design (`docs/10x_chemistry_reference.md`); its 8 sibling Multiome-GEX samples processed normally | protocol text `Single Cell Gene Expression Flex Fixed RNA Profiling (FRP)` | 18 | not measured |

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

**Batch 10, the same pathology at twice the rate** — 75 of 748 aligned samples (10.0%), 967
CPU-hours. The title vocabulary is wider, and CellPlex (`CMO`) is new:

| Accession(s) | What it actually is | Titled in GEO as | Samples |
|---|---|---|---|
| GSE219098 | Antibody capture | `Plasmablast v1D14 Rep1 - ADT`, … | 27 |
| GSE219153 | Barcode libraries | `BC1` … `BC12` | 12 |
| GSE218391 / GSE218392 | CellPlex multiplexing oligo | `CMO_NHCF_72h_IL1b`, … | 9 |
| GSE220730-220738 | Feature barcode, pooled tumours | `…_FB_pool_of_5_tumors_CFP_…` | 8 |
| GSE218483 | Hashing | `Monocytes 24h, HTO-derived cDNA` | 4 |
| GSE217689 / GSE217690 | Lineage/CRISPR barcode | `sc_rep_promoter_series_oBC_CS1_repB2`, `…mBC…` | 3 |
| GSE220085, GSE220651, GSE220468, GSE217976 | CMO / feature / ADT / HTO | `Pooled CMO library for 9 brain organoids`, `-ADT`, `hITM_HTO` | 7 |

For all but GSE218316 and the GSE218391/GSE218392 pair, the submission's GEX half processed
normally and only the feature half was wasted.

**Batch 11, the worst of the three — 144 of 844 aligned samples (17.1%)** — and the batch that
adds two classes the earlier two did not contain.

The familiar kind, 24 samples, caught by the same `exon_u` screen:

| Accession | What it actually is | Titled in GEO as | Samples | `all_u` | `exon_u` | Med_nFeature |
|---|---|---|---|---|---|---|
| GSE221776 | CITE-seq antibody capture | `PBT_batch1_library1_7donors_CITE`, … | 12 | 0.0001-0.0004 | < 0.00001 | 1 |
| GSE221853 | Cell hashing / multiplexing | `BC1` … `BC13` | 12 | 0.0017-0.0028 | 0.00006-0.0006 | 2-7 |

**New in batch 11: V(D)J libraries, which map well and look healthy.** 23 samples, all titled
`…_TCR`, aligned at 73-90% with up to 246 median features — so neither the `exon_u < 0.05` screen
nor a mapping-rate gate can see them. They are 10x 5' TCR-enriched amplicon libraries, and the
matrix says so directly: in `GSM6896067` 14 of the top 15 genes by total UMI are TRBV/TRAV
segments (`TRBV7-2` 4.4%, `TRBV20-1` 3.9%, `TRBV5-1` 3.4%, …), with `MALAT1` the only non-TCR
entry.

| Accession | What it actually is | Titled in GEO as | Samples | `all_u` | `exon_u` | Med_nFeature |
|---|---|---|---|---|---|---|
| GSE221776 | 10x 5' V(D)J (TCR) | `…_TCR` | 18 | 0.73-0.88 | 0.19-0.34 | 46-246 |
| GSE222011 | 10x 5' V(D)J (TCR) | `…-TCR`, `…_TCR` | 5 | 0.86-0.90 | 0.34-0.43 | 61-120 |

GSE221776 is the clean worked example of the whole class: 48 requested samples split exactly
18 `_GEX` / 18 `_TCR` / 12 `_CITE`, and only the 18 `_GEX` are gene expression. The suffix is the
only signal that separates them, and it is 100% accurate here.

**And GSE222431, which is not 10x at all** — 97 samples, 11.5% of the batch on its own. Bulk
immunoglobulin repertoire sequencing: 5'RACE PCR amplicons on an Illumina MiSeq, processed
upstream with MiGEC and IMGT HighV-QUEST, with no cell barcode anywhere in the design. Reads are
341 bp and 261 bp of variable length. It was accepted because chemistry inference has no upper
bound on the barcode read's length, so the 261 bp read was assigned the barcode role on a 28.6%
whitelist hit (a constant 5'RACE primer prefix) and reported `layout_confidence: "unique"`. All
97 published with exactly 1 median feature. Unlike the feature-barcode cases above, **the archive
metadata does separate this one**: `!Sample_library_selection = RACE` and
`!Sample_instrument_model = Illumina MiSeq` on all 150 samples. Tracked in
`issues/2026-09-20-batch11-barcode-read-length-unbounded.md`.

The rest of GSE215253 is the counter-example worth keeping in view: its 54 enrichment-PCR runs
*were* correctly rejected at `RENAME10XRUN`, and its 41 scRNA-seq samples all processed cleanly.
Only two of its feature libraries slipped through, because their biological read is 42 bp rather
than the 20 bp of the rest. No length threshold separates that from good data —
`GSM6660161` (56 bp), `GSM6668926`, `GSM6681083` and `GSM6681046` (55 bp) all map at 78-94%.

**Batches 14 and 15 — 96 of 814 (11.8%) and 68 of 780 (8.7%).** Both are dominated by the
familiar feature-barcode class, and batch 14 adds the token the scan had been missing.

| Accession | What it actually is | Titled in GEO as | Samples | `exon_u` | Med_nFeature |
|---|---|---|---|---|---|
| GSE229733 | Antibody-derived protein (CITE-seq) | `Patient <n> prot`, `Healthy donor number <n> HD<n> prot` | 61 | ~1e-6 | 1 |
| GSE230227 | 10x 5' V(D)J (TCR and BCR) | `PBMC, <id>, …, TCR` / `…, BCR` | 22 | 0.33-0.60 | 59-418 |
| GSE229626 | 6 x cell-surface protein, 6 x V(D)J | `Donor_<n>_CSP`, `Donor_<n>_VDJ` | 12 | 2e-7 / 0.18-0.23 | 1 / 180-564 |
| GSE233046 | 18 x ADT, 13 x TCR amplicon | `…, ADT-derived cDNA`, `…, TCR-derived cDNA` | 31 | 3e-7 / 0.29-0.50 | 1 / 146-664 |
| GSE231523 | 10x 5' V(D)J (BCR) | `Patient <n>, BCR-derived cDNA` | 14 | 0.10-0.46 | 25-138 |
| GSE233230 | 10x 5' V(D)J (BCR) | `C<nnn> B cell, young/elderly, VDJ (BCR)` | 11 | 0.45-0.58 | 72-422 |
| GSE231370 | Cell hashing | `Donor <n> Hashtag scRNAseq` | 3 | 4e-7 | 1 |
| GSE230563 | Cell hashing | `HTO_A`, `HTO_B` | 2 | 0.003 | 1 |
| GSE231319 / GSE231320 | CellPlex multiplexing | `control (Multiplexing Capture)`, `polCA (…)` | 2 | 5e-7 | 1 |
| GSE231688 / GSE231689 | Cell hashing | `…, HTO, scRNA-seq, 5'` | 2 | 4e-7 | 1 |
| GSE230629 | 10x V(D)J (TCR) | `…_VDJ_TCR_pooled`, `…_SP_TCR_pooled` | 2 | 0.52 / 3.5e-7 | 116 / 1 |
| GSE229402 | CITE-seq antibody capture | `HC_B_cells_ADT` | 1 | 1.2e-6 | 1 |
| GSE230612 | CellPlex multiplexing | `MCF7 - CMO Multiplexing Capture` | 1 | 1.3e-6 | 1 |

**`prot` was the gap, and it cost 61 samples in one dataset.** GSE229733 is 122 samples split
exactly in half: 61 titled `… geneexp`, which are healthy GEX, and 61 titled `… prot`, which are
the antibody arm. The `exon_u` screen caught all 61; the token rule below caught none, because
`prot` was not in it. The two-word phrase `Multiplexing Capture` was likewise unmatched.
**Both should be added to the pattern.** GSE229733 is otherwise the cleanest worked example in
this file — a 50/50 split where the suffix is the only signal and it is 100% accurate.

**Measurement, not the title, settled the marginal cases.** V(D)J-segment UMI fraction
(`TR[ABGD][VDJC]|IG[HKL][VDJC]` over the filtered matrix) against a GEX baseline of 0.28%
measured on `GSE232824`/`GSM7384003`: confirmed V(D)J libraries run 3.0-47.3%, and **13
title-scan hits across the two batches measured at or below baseline and are excluded from the
counts above** — GSE232822/GSE232824's 6 (`T cells, NY-ESO-1 TCR + library knockin`, where the
token names a knocked-in construct, not a library type, and which are healthy GEX at 2,313-3,634
median features), GSE230629/GSE230631's `GSM7230115` (`…_GEX_TCR_pooled`, the GEX arm of the
pair), and all 6 of GSE229279's, which are their own unresolved problem and have a row in the
Defects table above.

Two further batch-14 samples flagged by the `exon_u` screen are **not** non-GEX and are excluded:
GSE228957's `GSM7146182` and `GSM7146187`, both titled `…, single-nucleus`. Their matrices hold
real biology (`PAX5`, `PLCG2`, ribosomal genes) — they are shallow libraries where cell calling
over-called badly, 46,790 and 20,367 cells from ~400k and ~365k total UMIs, so ~8-18 UMI per
cell. Note this is *not* the batch-10 single-nucleus false positive pattern, which had
`full_u` 0.47; here `full_u` is 0.025, so the conjunction screen was right to flag them and only
inspection of the matrix separates the two.

**Batches 16 and 17 — 156 of 737 (21.2%) and 37 of 785 (4.7%).** Batch 16 is the worst
proportion measured so far, and the `exon_u` screen catches **3** of its 156.

| Accession | What it actually is | Titled in GEO as | Samples | `exon_u` | Med_nFeature |
|---|---|---|---|---|---|
| GSE234181 / GSE234185 | Lineage-barcode amplicon, paired with GEX of the same cells | `MCF7_AI1m1_BC.S31`, `T47D_AI2m1 _BC.S13`, … (GEX arm: `…_Cell.S<n>`) | 92 | 0.17-0.29 | 13-662 |
| GSE234069 | 10x 5' V(D)J (BCR and TCR), paired with GEX | `CR061_BCR_scRNA-seq`, `CR061_TCR_scRNA-seq` (GEX arm: `CR061_GEX_scRNA-seq`) | 32 | 0.17-0.63 | 9-827 |
| GSE233304 | 10x 5' V(D)J (TCR), paired with GEX | `Patient 15, Tumor CD45+, 5' TCR scV(D)Jseq` | 17 | 0.23-0.47 | 126-713 |
| GSE234778 | 8 × V(D)J, 4 × GFP amplicon | `D1_1_GFP+_VDJ`, `D1_1_GFP+_GFP` | 12 | 0.39-0.40 / 0.0009-0.053 | 119-250 / 1-8 |
| GSE233506, GSE233953 | 10x V(D)J (BCR) | `PBMC1_BCR`, `lymph node S144 scRNAseq BCR enrichted` | 3 | 0.19-0.35 | 44-135 |
| GSE237034 | Cell hashing | `Multipled-HashTag_Fibroblast Line A01 to A05, scRNAseq` | 8 | ~3e-7 | 1 |
| GSE236057 | CROP-seq guide amplicon | `NHA_CROPseq_Lib1_guidePCR` | 7 | 0.0002-0.0005 | 1-2 |
| GSE235663 / GSE235665 | 6 × CITE-seq, 6 × TCR | `014_10x_048_JiCh03_A_CITE`, `…_TCR` | 12 | ≤ 0.0004 / 0.10-0.35 | 1 / 28-83 |
| GSE235917 / GSE235920 | Not established; titled as V(D)J | `C1_mRNA_VDJ`, `Chir_Tils_mRNA_VDJ` | 4 | 0.32-0.35 | 26-39 |
| GSE235604 | 2 × CSP, 2 × V(D)J | `CMV CSP, 10 donors, scRNAseq`, `4virus VDJ, …` | 4 | ~1e-7 / 0.20-0.24 | 1 / 99-372 |
| GSE235072 | Cell hashing | `…, hashtag oligos library 1` | 2 | ~1e-6 | 1 |

**The lineage-barcode class is new, and neither screen sees it.** GSE234181 is a 184-sample
series split exactly in half between `…_Cell.S<n>` (healthy GEX) and `…_BC.S<n>` (the barcode
amplicon). The amplicon maps at `exon_u` 0.17-0.29, far above the 0.05 screen. Only 58 of the 92
fall under 100 median features, and no title token matches `_BC.S13`. `GSM7453678` holds 140,768
UMIs across 10,410 cells, about 14 per cell. Add `_BC.S<n>` to the pattern.

**Every V(D)J library in batch 16 has a GEX sibling of the same cells, and the sibling is the
best evidence.** GSE233304's `5' TCR scV(D)Jseq` libraries hold 0.15-1.1 M UMIs where their
`5' GEX` siblings hold 8-45 M. **The V(D)J-segment UMI fraction does not separate them here**,
unlike in batches 14-15: the GEX siblings measure 0.0-5.8% and the TCR libraries 0.4-8.7%, so
the ~3% decision line above is not general. The title plus the paired sibling is what settled
these.

**Batches 18 and 19: 49 of 895 (5.5%) and 67 of 767 (8.7%).** 5 samples are in both
batches, reached through two dataset accessions, so 111 are distinct. The `exon_u` screen
catches all 64 feature-barcode libraries and none of the 47 V(D)J libraries. Every V(D)J library
here has a GEX sibling of the same cells, and every one holds 10-50× fewer UMIs than that
sibling (0.1-6.7 M against 3.8-85 M). The V(D)J-segment UMI fraction has a median of 20%,
against 0% in the siblings, but 12 of the 47 are under 10%. The sibling comparison settled
these, as in batch 16.

| Accession | What it actually is | Titled in GEO as | Samples | `exon_u` | Med_nFeature |
|---|---|---|---|---|---|
| GSE240252 | 10x V(D)J, paired with GEX | `BL-101_VDJ` (GEX arm: `BL-101_GEX`) | 15 | 0.14-0.48 | 83-374 |
| GSE243002 | 10x V(D)J (TCR and BCR), paired with GEX | `IgAV_2_BCR`, `control_1_TCR` (GEX arm: `…_rnaseq`) | 12 | 0.13-0.42 | 75-353 |
| GSE241842 | 10x V(D)J (TCR), paired with GEX | `CML patient 5 TCR` (GEX arm: `… CD3`) | 8 | 0.19-0.26 | 35-387 |
| GSE240865 | 6 × V(D)J, 6 × ADT | `Memory CD4+ T cells …, control - VDJ`, `… - ADT` | 12 | 0.38-0.53 / ~1e-7 | 54-218 / 1 |
| GSE238120 | 4 × BCR, 4 × feature barcode | `B cell, DT-sorted, BCR`, `… Feature Barcode` | 8 | 0.49-0.54 / ~1e-6 | 90-145 / 1 |
| GSE242775 | 10x V(D)J (BCR), paired with GEX | `post-vaccination Delta vdj library` | 2 | 0.55-0.57 | 375-586 |
| GSE238189 | Feature-barcode capture | `GE_p38i_rpl2__FBC` | 8 | 0.0009-0.008 | 1-25 |
| GSE243006 / GSE243005 | CITE-seq ADT and HTO | `scRNAseq_CITEseq_CD34+cells_CTRL1_ADT` | 11 | 0.00008-0.0034 | 1-3 |
| GSE243756 / GSE243765 | CROP-seq guide libraries | `grna library 1, scRNA-seq` | 12 | 0.0003-0.0006 | 2 |
| GSE243443 / GSE243445 | Cell hashing | `Control colonoids from C9, HTO` | 7 | 0.0007-0.003 | 1-18 |
| GSE243349 | Feature-barcode capture of hashtags | `Spleen, combined and hashtagged, #12, #17, #18 FBC` | 5 | ~1e-6 | 1 |
| GSE241683 / GSE241838 / GSE241839 / GSE241882 | Hashtag, MULTI-seq and CROP-seq guide libraries | `Hashtag T cell Cas9 genome edited`, `…, CROP-seq gRNA` | 4 | 1e-5-0.012 | 5-37 |
| GSE243501 | CellPlex multiplexing | `Multiplexed_CMO: [468, CSF] + …` | 3 | ~6e-7 | 1 |
| GSE242425 / GSE242426 | Feature barcode | `D20organoids-FB` | 2 | 0.00007-0.0001 | 1-3 |
| GSE241739 / GSE242039 | CellPlex multiplexing | `Subject_5_Blood = 303, Subject_5_Tissue = 304 (CMO)` | 1 | 1.3e-6 | 1 |
| GSE240441 | Cell hashing | `Tonsils and PBMC NK cells, HTO` | 1 | 2e-7 | 1 |

**Two title forms the current token list misses.** `(CMO)` is delimited by parentheses, which
are not in the delimiter class, so GSE241739's CellPlex library passed the title scan. It was
caught only by the `exon_u` screen. `FBC` and `-FB` (feature-barcode capture) are not tokens
at all. Add `()` to the delimiter class and `FBC`/`FB`/`MULTI-seq` to the token list.

Three GEX-titled samples also sit under the screen with no explanation found: GSE239750's
`tonsil tissue-SE-GEX` (908 k reads, 82 cells, `exon_u` 0.0006 but `full_u` 0.84),
GSE239889's `248N` (562 M reads, 1.5% in cells) and GSE240766's `H1 PC C` (1.0% in cells).
They are held as suspects, not counted above.

Excluded as false positives of the widened scan: GSE237239's `mRNA_Hashtag_R342_DMSO` (the GEX
arm, 2,392 median features) and GSE236107's `Mouse_T1_GFP` etc. (GEX from a GFP sort). A bare
`GFP` token is unsafe; `_GFP$` is not.

**Batch 22: 190 of 1,838 (10.3%)**, 116 of them V(D)J. Every V(D)J library has a GEX sibling in
the batch; the V(D)J libraries hold 0.07-2.0 M UMIs and 3.2-99.7% of them in TR/IG segments,
against 28-140 M and 0.1-2.3% in six measured siblings.

| Accession | What it actually is | Titled in GEO as | Samples | `exon_u` | Med_nFeature |
|---|---|---|---|---|---|
| GSE249131 | 10x V(D)J (BCR and TCR), paired with GEX | `Healthy1_PBMC_BCR`, `VEXAS9_PBMC_TCR` (GEX arm: `…_GEX`) | 44 | 0.21-0.56 | 7-103 |
| GSE248556 | 10x V(D)J (TCR and BCR), paired with GEX | `A_LC11_TCR`, `A_LC11_BCR` (GEX arm: `A_LC11`) | 30 | 0.20-0.62 | 82-467 |
| GSE253352 | 10x V(D)J (TCR), paired with GEX | `Biopsy IL15.CAR TCR Patient 2` (GEX arm: `… 15.CAR mRNA Patient 2`) | 22 | 0.25-0.53 | 24-474 |
| GSE252331 | Antibody-derived tags | `Sepsis 5, CCI, Day 14, ADT` | 17 | ~2e-7-8e-7 | 1 |
| GSE250242 / GSE250243 | 9 × V(D)J, 9 × CSP | `PSA1 VDJ`, `PSA1 CSP` | 18 | 0.24-0.49 / ~3e-8-3e-5 | 28-148 / 1 |
| GSE250235 | CellPlex multiplexing | `Day 0 CMO` | 10 | ~2e-7-3e-6 | 1 |
| GSE250378 | CROP-seq guide libraries | `Day 0, Replicate 1, CROP` | 9 | 0.004-0.006 | 5-8 |
| GSE252416 / GSE252642 / GSE252687 | 10x V(D)J (BCR), paired with full transcriptome | `lymph node … scRNAseq BCR enriched` | 7 | 0.23-0.41 | 9-229 |
| GSE251912 | Cell hashing | `Islet_150 (SAMN21845647), HTO` | 5 | 6e-5-0.0022 | 1-45 |
| GSE252724 | Cell-surface protein | `SJ donor 1 and SJ donor 3 Biopsy (…) - CSP` | 5 | ~2e-7-4e-7 | 1 |
| GSE249006, GSE249597, GSE248788, GSE158702, GSE248951 | Hashtag / HTO, CROP-seq gRNA | `antibody hashtag …`, `Library_1, HTOs`, `Library_1, gRNA`, `Hash tag oligo library 1`, `HTO reads` | 13 | ≤ 0.007 | 1-16 |
| GSE248788, GSE253828 | 10x V(D)J (TCR) | `scRNAseq, TCR library 1`, `SP019-H,scTCR-seq` | 4 | 0.29-0.48 | 21-57 |
| GSE248590, GSE252830 | Barcode libraries | `Mock KO CART19 BC`, `… T cell scRNAseq Barcode` | 4 | ~4e-7-1e-6 | 1 |
| GSE153931 | CITE-seq antibody capture | `CoVi17_8_S_24_9P_CITE` | 1 | 0 | no filtered matrix |
| GSE248549 | Not established | `TLB, scRNAseq` | 1 | 0.021 | 8 |

**The title scan's largest false positive: GSE249313, 44 samples.** Every GEX sample is titled
`scRNA-seq and TCR profiling of Homo Spaiens: <patient> <day> …`, and the `TCR` token matches
the study description. They are healthy GEX at 1,457-2,380 median features. The scan also
missed 20 of the 190: `CROP` (9), `Hash tag oligo` (3), `HTOs` (2), bare `BC` and `Barcode`
(4), `scTCR-seq` (1) and the untitled `TLB` library. The `exon_u` screen caught all of those
except `scTCR-seq`.

Tracked, with the measured gate that would catch them, in
`issues/2026-09-18-batch9-no-mapping-rate-gate.md`; the V(D)J class, which that gate cannot
catch, in `issues/2026-09-20-batch11-vdj-libraries-published-as-gex.md`. Nothing in this section
is worth reporting to GEO — every submission here is accurate about what it is.

A delimited-token rule on `!Sample_title` (`TCR`, `BCR`, `VDJ`, `CITE`, `HTO`, `ADT`, `CSP`,
`CMO`, `hashtag`, `hashing`, `gRNA`, `sgRNA`, `prot`, the phrase `Multiplexing Capture`, plus
bare `BC<n>`) flags 47 of batch 11's 844
samples with **zero false positives** — every feature-barcode and V(D)J sample above. Token
delimiting is what makes it safe: a substring match instead pulls in `Ascites_IG_1` (a specimen
code), `scRNA CRISPRa single perturbation S1` (the GEX arm of a Perturb-seq study, 3503 median
features) and `spike-specific B cells, D14, scRNA`, all of which are healthy GEX.

**On batches 14 and 15 the same rule is no longer false-positive free**, and the exceptions are
instructive: it flagged 13 samples that measurement cleared (above). A `TCR` token can name a
knocked-in construct or the GEX arm of a paired submission as easily as a V(D)J library, so
treat a title hit as a candidate and confirm it against the matrix — the V(D)J-segment UMI
fraction separates them cleanly, and a GEX baseline near 0.28% makes anything above ~3%
decisive.

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
