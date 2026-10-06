# Baseline distributions and open issues

Distilled from `docs/agent_debug.md` §7 and §13, plus what the September 2026 collections added.

## §baseline — the August 2026 run (1814 failures)

A reference distribution. If a new run's profile differs sharply from this, that difference is
itself the finding.

Stage: `RENAME10XRUN` 1738 · `RENAME10XSAMPLE` 39 · `STARSOLO10X` 23 · `WGET10X` 14.
Totals: 1800 classified + 14 with no log block = 1814.

| Cause | Runs | Verdict |
|---|---|---|
| `chem-no-whitelist-hit` | 1014 | 880 genuinely not 10x (GSE109816); the rest mixed |
| `bc-variable-length` | 401 | pipeline too strict — chemistry called at 96%+ |
| `no-single-gex-r2` | 159 | correct: feature-barcode (ADT/HTO) libraries, R2 ≈ 25 bp |
| `single-fastq-from-sra` | 89 | `SRA2FASTQ` emitted one mate |
| `read-counts-differ` | 48 | `SRA2FASTQ` partial dump |
| `r1-len-differs-runs` | 39 | missing chemistry JSON wiring |
| `chem-no-whitelist-rand` | 13 | STARsolo re-detect disagreeing with run-level inference |
| `chem-ambiguous-layout` | 13 | too strict |
| `truncated-gzip` | 9 | download integrity |
| `star-sj-buffer` | 6 | `--limitOutSJcollapsed`; 5 at exit 104, 1 at 139 |
| `corrupt-fastq` | 5 | malformed record mid-file |
| `UNCLASSIFIED` | 4 | 3 = LSF `TERM_MEMLIMIT` at exit 130; 1 = a `starsolo_10x_auto.sh` check |
| no log block | 14 | `WGET10X`: 13 × HTTP 403 (transient), 1 × empty URL |

Roughly 1040 correct rejections, 774 avoidable.

## §later — what the later batches looked like

Failure volume dropped by two orders of magnitude once the guards landed. The remaining
failures are a different population, so do not expect the August shape.

| Batch | Run | Failed tasks | Permanent | Notes |
|---|---|---|---|---|
| 3 | `big_keller` | 111 | 92 | 41 UNCLASSIFIED — that manifest predates the LSF-kill line, so exit 130/139 carry no evidence |
| 5 | `tender_brattain` | 48 | 48 | 32 `sra-zero-length-read` in one dataset (GSE200629 = 24 of 48) |
| 6 | `spontaneous_ampere` | 16 | 7 | 9 self-healed on retry; 3 × `chem-no-whitelist-rand`, 2 × `star-segfault`, 1 × `meta-no-run-id`, 1 × `sra-too-few-fastq` |
| 7 | `elegant_lamarr` | 10 | 9 | 1 self-healed (`lsf-memlimit`); 8 × `chem-ambiguous-layout` (one dataset, GSE212038 — truncated-UMI allowlist gap, see `§open`), 1 × `chem-no-whitelist-hit` (genuinely near-empty source FASTQ, confirmed against the raw download). First batch with zero datasets lost. |
| 8 | `cheeky_cuvier` | 53 | 11 | 42 self-healed (35 × `wget-exhausted`, 6 × `truncated-gzip`, 1 × `lsf-memlimit`); **10 of the 11 permanent are one pathology in one dataset pair** — GSE212964/GSE212965 split-mate runs, 5 samples, see `§open` — plus 1 near-empty SRA object (GSE213370). 6 of 787 unique samples lost. Also the first batch where duplicate dataset rows were quantified: 200 of 981 alignment tasks redundant (20.4%). |
| 9 | `gigantic_hypatia` | 108 | 102 | 6 self-healed (3 × `truncated-gzip`). **The failure list is the small half of this batch.** 68 of the 102 are correct rejections (54 ADT/CRISPR runs in GSE215253 alone); the two real defects are 9 samples lost to a STAR segfault in GSE214532 (`§open`) and 2 to a `umi_len` disagreement in GSE215362. Separately, and invisible in the failure list: **38 aligned samples are not GEX libraries at all** and were published as near-empty matrices, and 133 of 917 alignments (14.5%, 2131 CPU-h) were duplicates. 55 of 839 samples lost. |
| 10 | `trusting_carson` | 16 | **13** | A quiet failure list and the batch's worst outputs. 3 self-healed; the 13 span 4 datasets, of which 7 are GSE217689's barcode libraries. **The findings are all outside the failure list**: 75 of 748 aligned samples (10.0%, 967 CPU-h) are hashing/ADT/CMO/feature libraries published as near-empty matrices, 6 samples were dropped inside `FETCH10XMETA` at exit 0 with their reads pooled into a surviving sample's matrix (`§open`), and 225 of 973 alignments (23.1%, 4056 CPU-h) were duplicates. 14 of 762 samples lost. |

| 11 | `irreverent_elion` | 22 | **8** | The cleanest failure list of any batch and the worst output quality measured. 14 self-healed (11 × `lsf-memlimit`, 3 × `truncated-gzip`); the 8 permanent are 6 × GSE220774 `concatenated-read` (a *third* cause of `no-single-gex-r2`, see `§open`), 1 LSF memory kill, 1 threshold near-miss. **Everything that matters is outside the failure list**: 144 of 844 aligned samples (17.1%) are not GEX — 97 of them a bulk Ig amplicon dataset admitted by an unbounded barcode read length, 23 V(D)J libraries that map at 73-90% and defeat the §4 screen entirely, 24 CITE/hashing. 148 of 992 alignments (14.9%, ~1421 CPU-h) were duplicates. `FETCH10XMETA` emitted 852/852 — batch 10's silent drop did **not** recur. 8 of 852 samples lost. |
| 18 | `compassionate_shaw` | 27 | **16** | Ran concurrently with batch 19; the two share datasets under different accessions. 11 self-healed (3 × `lsf-memlimit`, 3 `WGET10X` killed externally, 1 truncated gzip, 2 `SRA2FASTQ`, 1 × STAR segfault, 1 argparse error on a mate-less run). 12 of the 16 are index-only runs registered separately from their libraries (GSE241739). Harmless: every sample aligned from its real runs. 3 are SRA runs with 0 spots (GSE241292). 1 is `GSM7744878`: 8.06 × 10⁹ reads, segfaulting in STARsolo's Velocyto step after writing a good Gene matrix (`§open`). **The batch-10 `FETCH10XMETA` shared-BioSample drop recurred**: GSE241292, 5 samples lost at exit 0 and 5 published matrices pooling two libraries. 49 of 895 aligned (5.5%) not GEX, 25 of them V(D)J the `exon_u` screen cannot see. 93 of 988 alignments (9.4%, 903 CPU-h) duplicates. 6 of 901 samples lost. |
| 19 | `focused_goldberg` | 24 | **16** | 8 self-healed (4 × truncated gzip, 2 argparse errors on mate-less runs, 1 × `lsf-memlimit`, 1 `SRA2FASTQ` at 140). 12 of the 16 are the **same** GSE241739 index-only runs, reached as GSE242039, and 1 is the **same** `GSM7744878` Velocyto segfault, reached as GSE241882. **New: 3 samples lost because ENA's two-file FASTQ omitted the 28 bp barcode read that SRA holds** (GSE241998/GSE241999, `§open`). A retrospective scan found batch 14's GSE228428 (5 samples) was the same defect, misfiled as a submitter one. 67 of 767 aligned (8.7%) not GEX, 22 of them V(D)J. **213 of 980 alignments (21.7%, 1,190 CPU-h) duplicates**, plus 16 samples aligned in both batches. 4 of 771 samples lost. |
| 20 | `sleepy_curie` | 104 | **94** | Ran concurrently with batch 21. 10 self-healed (6 × `lsf-memlimit`, 4 × truncated gzip). **52 of the 94 are the strand-test segfault** in one snRNA-seq dataset (GSE245310), at margins of 9-23 points, so the drafted margin fix would not have saved them. The test reads GeneFull; exonic Gene said Forward in all 52 (`§open`). 24 are a four-way split-run deposit (GSE245998, 6 samples) caught at `SRA2FASTQ` and verdicted correct-rejection by `triage.py`, wrongly. 18 are a correct rejection (GSE245175, SORT-seq). 64 of 577 aligned (11.1%) not GEX, 51 of them V(D)J. 131 of 708 alignments (18.5%, 1,997 CPU-h) duplicates. 67 of 644 samples lost. |
| 21 | `sick_galileo` | 30 | **27** | A short failure list and the batch's biggest loss outside it: **`FETCH10XMETA` dropped 159 of 873 samples at exit 0** (GSE246613 147, GSE247111 12, shared BioSample) and published 112 pooled matrices, most under the TCR/CellPlex GSM's name. 3 self-healed (truncated gzip). The 27: 18 split-mate (GSE247205, 9 samples), 8 truncated-UMI `len=27` (GSE247827, whole dataset; `10,11` dry-run accepts it), 1 0-spot SRA run. 123 of 701 aligned (17.5%) not GEX, 108 of them V(D)J (GSE247531 alone 89). 114 of 815 alignments (14.0%, 1,245 CPU-h) duplicates; 24 samples also aligned in batch 20 (GSE245187 = GSE246960). 172 samples lost, all recoverable. Work dirs for batches 6-19 found deleted this session. |
| 22 | `gloomy_church` | 73 | **29** | **The first batch from `batches_deduplicated/`, and 0 duplicate alignments** (1,838 QC rows = 1,838 samples). 44 self-healed: **35 × STARsolo `lsf-memlimit`** at 128 GB, every one a sample of ≥ 1.35 × 10⁹ reads, 2,686 core-h (7.6% of the stage), plus 9 truncated gzip. Of the 29, 19 are harmless (18 index-only runs in GSE248788/GSE249894/GSE250444, one 2,686-spot run), which `triage.py` leaves as "investigate". **New: SRA runs whose read layout changes partway through** (lanes with 3 vs 2 reads per spot, or R1/R2 loaded as consecutive single-end spots) fail `SRA2FASTQ`'s mate-count guard with the data complete — 3 samples (GSE249159 ×2, GSE248489 ×1, `§open`). **New: ENA FASTQs corrupt at source** — 4 files of 2 GSE249894 runs fail `pigz -t` while matching ENA's own `fastq_md5` (`§open`). The shared-BioSample drop recurred (GSE254170, 1 lost, 1 pooled). A mate-less run cost data for the first time (`GSM7996301`, 7.7% of reads), because deduplicated input no longer supplies a rescuing duplicate. **190 of 1,838 aligned (10.3%) not GEX**, 116 V(D)J; the title scan flagged 44 healthy GEX (GSE249313, `scRNA-seq and TCR profiling …`) and missed 20. 4 of 1,842 samples lost, 3 published partial, all recoverable. |

Batch 6 is the useful worked example: `.nextflow/history` records it `OK`, and it still had 7
permanent failures. `errorStrategy 'ignore'` means a green run tells you nothing.

Batch 11 is the other one worth internalising: 22 failed tasks, 8 permanent, no dataset lost —
and it published more bad data than any batch before it. A short failure list is not good news,
it is a reason to run §4.

## §open — found but not fixed

* **A sample whose runs mostly failed is published anyway, from whatever survived.**
  `RENAME10XRUN` runs per run; a refused run is dropped and the sample continues with the rest,
  so `RENAME10XSAMPLE` merges a short group and a matrix is published under the sample's name
  with nothing recording that it is partial. Batch 12-13: **5 samples**, and **three of the five
  look entirely healthy** — `GSM7056649` was built from **1 of its 17 runs** and came out at 1262
  median features, 18569 cells, 39% exonic. Only `GSM7056650` (7 median features) is visible to
  the §4 near-empty screen. Same class as batch 10's `FETCH10XMETA` drop but one stage later, and
  unlike batch 10 the evidence is already in the run's own artefacts: the expected run count is
  the `GroupKey` size and is never compared against what arrived. Drafted in
  `issues/2026-09-23-batch12-13-partial-run-merge-silent.md` — no diff, the choice between
  recording, warning and gating is a policy decision. **The four-line join that detects it belongs
  in `verify.md §qc` regardless of which option is taken.**
* **`INDEX_LENGTHS` is an exact-length allowlist and has no 12.** In
  `infer_10x_run_recommended.py:157`, `assign_index_roles` refuses a run whose leftover files are
  not of an allowed index length. GSE225807's four-file SRA dumps carry a constant 12 bp i7
  beside an 8 bp i5, a 16+10 barcode read matching the 3'v2 whitelist at 94.8%, and 258 bp of
  cDNA — valid 10x 3' GEX, rejected on 20 runs, **the largest single cause in batch 12-13** (20 of
  65 permanent). The guard also names the whole leftover list rather than the offenders, so the
  message blamed the acceptable 8 bp file too. Diff attached in
  `issues/2026-09-23-batch12-13-index-length-12-rejects-run.diff` (adds `12`, narrows the
  message); the `12` itself is a judgment call, and the allowlist design will need widening again.
  `infer_10x_run.py:155` carries the same constant and is not the wired script — settle whether it
  is dead code, as the two `sra2fastq` modules were.
* **Samples with zero gene-mapped reads are published with no filtered matrix, and the QC table
  misreports them** (GitHub #17; corrected 2026-10-02). Six batch 12-13 samples (GSE224986:
  `GSM7036555`-`559` and `561`; `560` is a normal row) show `WL = Undef`, **6,794,880 cells**,
  1 median feature, `exon_u` 0. This was first read as "aligned with no whitelist", and that
  was wrong. The driver log shows all six got `--wl gex_3pv3_family`. The values are what
  `starsolo qc` prints from its plate-based branch when `output/Gene/filtered/` is missing:
  `Undef` is hard-coded and 6,794,880 is the raw v3 barcode list, not a cell count.

* **Chemistry inference has no upper bound on the barcode read length.**
  `gex_observed_umi_len()` in `infer_10x_run_recommended.py` accepts any read where
  `observed = len - cb_len >= chem.umi_len`, so a 261 bp bulk amplicon read is accepted as a
  28 bp GEM-X v4 barcode read. Every length option the script exposes is a minimum; there is no
  maximum, and no way to impose one from `ext.args`. The ATAC path in the same file *does* have
  the bound, as `info.common_length in ATAC_BARCODE_LENGTHS`. `require_constant()` is entered
  (the read is only 47.5% constant against a 0.98 threshold) and then accepts it as
  "quality-trimmed", because its tests are lower bounds too. Cost: batch 11's GSE222431 — bulk
  5'RACE Ig repertoire data on a MiSeq, not 10x in any sense — passed on all 97 runs with
  `layout_confidence: "unique"` and published 97 matrices of 1 median feature. Across those 97
  runs the length assigned the barcode role was 260-271 bp. Fix is a few lines but sits in a
  registry module with a `.module-info` checksum, so no diff is attached; drafted in
  `issues/2026-09-20-batch11-barcode-read-length-unbounded.md`.
* **V(D)J libraries map well and cannot be caught by any mapping-rate gate.** Batch 11 published
  23 TCR-enriched 5' amplicon libraries as GEX (GSE221776 ×18, GSE222011 ×5). They align at
  73-90% with 46-246 median features and look entirely healthy; `GSM6896067`'s matrix has 14 of
  its top 15 genes as TRBV/TRAV segments. The §4 screen catches 6 of 23 and the `max(fwd,rev)
  < 5%` gate proposed for the feature-barcode case catches none — both key on the library
  mapping to *nothing*. The only signal is the GEO title. A delimited-token rule on
  `!Sample_title` flags 47 of batch 11's 844 samples (all 23 V(D)J plus the 24 CITE/hashing) with
  **zero false positives**; a substring match instead of a token match gives 13. Drafted in
  `issues/2026-09-20-batch11-vdj-libraries-published-as-gex.md` — no diff, rejecting on
  submitter free text is a design decision.

* ~~**`errorStrategy` precedence bug.**~~ **Closed — not a bug.** Retrying every task once
  whatever its exit code is intentional, and the `task.attempt == 1` clause is load-bearing:
  batch 6 lost no sample to a transient exit-1 failure because of it. See `classify.md §exits`
  for the evidence table. Two things remain true and are not defects: `attempt` is not
  diagnostic, and a deterministic failure costs one extra attempt — the accepted price of not
  discarding transient ones.
* **No failure manifest inside the pipeline.** `bin/collect_run_logs.sh` reconstructs one after
  the fact, but a `workflow.onComplete` hook writing run/sample/dataset/process/exit/first-error
  would remove the whole collection step. A previous attempt at this was reverted (commit
  `028aaf3`) for causing silent failures.
* **`--allow-gex-truncated-umi-lengths` default (`10`) is missing `11`.** In
  `infer_10x_run_recommended.py`, this flag is the guard's own named mechanism for a v3 GEX read
  that is short on UMI bases, but its default only accepts an observed UMI of 10. Batch 7's
  GSE212038 (8 runs, `SRR21201139/156/157/158/159/164/165/166`) has R1 uniformly 27 bp — 16 bp CB
  + 11 bp UMI — which falls through the allowlist and gets the generic "no supported 10x run
  layout matched" rejection instead of the truncated-UMI acceptance path. Confirmed by re-running
  the script against the persisted work dir with `--allow-gex-truncated-umi-lengths 10,11`: it
  then accepts the run with a truncated-UMI warning, no other change. Fix is one value, in one
  place — `ext.args` for `RENAME10XRUN` in `configs/rename10xrun.config`. Not yet applied; see
  the [Batch 7 Post-Mortem](../../../../docs/post-mortems.md) for the full evidence. **Still unapplied at
  2026-10-02, and batch 21 lost a whole dataset to it** (GSE247827, 4 samples, 8 runs at
  `len=27`, 95.7-96.4% v3); the `10,11` re-run against its work dir accepts it. 8 samples to date.
* **Split-mate runs are rejected per run, upstream of the merge that would fix them.** When a
  submitter registers R1 and R2 of one 10x library as two SRA run accessions, `RENAME10XRUN`
  infers layout per run and rejects both halves — the barcode half for having no biological
  read, the cDNA half for having no whitelist hit — although `links.tsv` already names the same
  GSM on both rows and `RENAME10XSAMPLE` exists downstream to merge them. Cost so far: 5 samples
  in batch 8 (GSE212964/GSE212965, measured) and 1 in batch 5 (GSE202476, log-text signature
  only — its work dirs are gone). The detection signature is a sample whose runs fail with
  *complementary* errors. Fix needs a design decision, not a threshold; drafted in
  `issues/2026-09-16-batch8-split-mate-runs.md`. **Batches 20-21 added 15 samples**:
  GSE247205 (9, the same two-way shape) and GSE245998 (6, a **four-way** split, I1/I2/R1/R2
  as four SRA runs). The four-way split fails at `SRA2FASTQ`, a stage before any fix proposed in
  the issue. Better signature, measured: **runs of one SRX with identical spot counts**, read
  from `sra.tsv`. Across both batches' 790 multi-run SRXs it matched exactly the 15 split ones.
* **Datasets with identical sample lists are processed twice.** GEO SuperSeries/SubSeries pairs
  whose `sample_id` field is byte-identical become two dataset units, so every shared sample is
  downloaded, renamed and aligned twice. Batch 8: 200 of 981 alignment tasks redundant
  (20.4% of the stage, 26% more than the 781 distinct samples needed), plus 545 redundant
  `WGET10X` and 56 redundant `SRA2FASTQ`. Batch 9: 21 duplicated pairs of 118 dataset rows, 133
  of 917 alignments (14.5%) and 2131 of 14 844 CPU-hours. Project-wide, 445 of 4602 dataset rows
  (9.7%) are exact duplicates. Not a correctness bug — the outputs are right — but the largest
  single waste measured so far, and it is not improving: batch 10 is 225 of 973 (23.1%) and 4056
  of 16 097 CPU-hours, the worst of the four measured; batch 11 is 148 of 992 (14.9%) and ~1421
  of 16 976 CPU-hours. It multiplies failures too: GSE214532's
  segfault burned 54 attempts rather than 27 because the dataset ran twice. It also means a per-sample tally taken
  off `mapping_qc_stats.tsv` without de-duplicating is over by 17%. Drafted in
  `issues/2026-09-16-batch8-duplicate-dataset-rows.md`.
* **`triage.py` verdicts `no-single-gex-r2` as a correct rejection.** The note asserts an
  ADT/HTO library, but the split-mate barcode half emits the same error text, and in batch 8 all
  5 rows were the latter — so the triage output read 54.5% correct-rejection for a set that was
  really 1 correct rejection and 10 tasks of one defect. Batch 5's 7 rows were misverdicted the
  same way (4 of them the GSE202476 sample). **Batch 11 found a third cause**: GSE220774's 6 rows
  are a *concatenated* read — CB+UMI+polyT+cDNA in one 150 bp record, with the 8 bp index as
  `_2` — and each affected sample has exactly one run, so the complementary-error signature
  proposed as the better fix would mislabel them as ADT/HTO. Any refinement now has to separate
  three cases, not two. One-line fix, dry-run confirmed, drafted in
  `issues/2026-09-16-batch8-triage-verdict-split-mate.md`.
* **`collect_run_logs.sh` writes a ragged trace when `nextflow log` warns on stdout.** Nextflow
  puts its FTP-proxy warnings on stdout, which is redirected straight into `runlogs<N>.tsv`; one
  stray line and duckdb's dialect sniffer refuses the whole file, with an error that mentions
  only CSV dialects. It bites from an agent session, whose sandbox sets `FTP_PROXY` to a
  `socks5h://` URL that Nextflow will not parse — i.e. exactly how this skill is normally run.
  All surviving traces (batches 6-8) are clean at 36 fields. Drafted, with diff, in
  `issues/2026-09-16-batch8-collect-trace-ragged.md`.
* **`FETCH10XMETA` can silently drop samples that share a run set — and pool their reads.** When
  several GSMs map to the same SRA runs, the metadata stage emits only the last of them, exits 0,
  and says nothing; the survivor is then aligned against *all* the group's runs. Batch 10's
  GSE218936: 8 requested, 2 emitted, 6 lost, and the two surviving matrices each pool four
  experimental conditions (1G/µG × stim/unstim) under one sample's name at 96% mapping and 1865
  median features — outputs that look perfect and are wrong. Cell count is not a detector (32709
  against a batch median of 7471 is only ~p97). The check that finds it is requested-vs-emitted,
  now in `SKILL.md §4`. Batch 9 is clean on it, and **batch 11 is too** — 852 of 852 unique
  requested samples emitted, so this did not recur; batches 1-8 cannot be checked.
  **It recurred in batch 18** (GSE241292, 2026-09-30): each "main sample"/"sub-sample" GSM
  pair shares one BioSample, so 5 main samples were lost at exit 0 and each of 5 published
  sub-sample matrices pools both libraries (`GSM7720741`: 161 M own reads + 691 M from the
  dropped GSM). The SRX-lookup diff drafted after batch 16 was still unapplied. Same family as
  GSE247111's `shared-biosample` in `docs/archive-pathologies.md`, different consequence. Drafted
  in `issues/2026-09-19-batch10-metadata-drops-samples-silently.md` — no diff, it needs a design
  decision about what a shared run set means. **Batch 21 is the worst yet** (2026-10-02):
  GSE246613 emitted 104 of 251 and GSE247111 10 of 22, so **159 samples (18.2% of the batch)
  were lost at exit 0**. 112 published matrices pool a GEX library into a TCR- or
  CellPlex-named GSM. The SRX lookup in the (still unapplied) diff resolves 273/273 to one run
  each. **Batch 22 again** (GSE254170): `GSM8427071` (uninduced control) dropped, its run
  pooled into `GSM8427072` (dox-induced knockdown), both on BioSample `SAMN42510873` with their
  own SRX. Running total: 186 lost, 125 pooled. GitHub #6.
* **STARsolo discards a finished Gene matrix when Velocyto segfaults.** `platform_10x.sh`
  hard-codes `--soloFeatures Gene GeneFull Velocyto`. On `GSM7744878` (8.06 × 10⁹ reads, 12
  runs, GSE241683/GSE241882) STAR writes `Gene` and `GeneFull` (106,104 cells, 3,876 median
  genes) and then dies at `Velocyto counting: allocated arrays`, exit 139. That happened on
  every attempt that got that far: 6 attempts across batches 18 and 19, ~1,360 core-hours. It
  cannot be worked around from `ext.args`: STAR rejects the repeated `--soloFeatures` with
  `FATAL INPUT ERROR: duplicate parameter` (measured). Needs a container flag. Drafted in
  `issues/2026-09-30-batch18-19-velocyto-segfault-discards-good-matrix.md`.
* **ENA's two-file FASTQ can omit the barcode read.** For a run SRA stores as 3-4 reads per
  spot, ENA exports two files, and it can drop the 28 bp barcode read (`SRR25819675`: SRA
  holds 91, 28, 8; ENA serves 91 + 8). `parse_metadata.sh` prefers `ENAFQ` whenever a `_1`/`_2`
  pair exists, so the run is rejected at `RENAME10XRUN` with `no whitelist hits` and looks
  like a correct rejection. Batch 19: 3 samples (GSE241998/GSE241999). Batch 14's GSE228428
  (5 samples) is the same defect and had been filed as a submitter defect. **Detection
  signature**: an `ENAFQ` run with one cDNA-length file and only index-length (≤ 10 bp)
  others. Check SRA's `sra-db-be/run_new?acc=<run>` read structure before calling it
  unrecoverable. Drafted in `issues/2026-09-30-batch19-ena-fastq-omits-barcode-read.md`.
* **A run with a permanently failed mate still reaches `RENAME10XRUN`, and fails with an
  argparse error.** When one `WGET10X` of a pair is ignored after 5 attempts, the surviving
  file goes on alone and `infer_10x_run_recommended.py` exits 2 with
  `This run-level script expects 2-4 FASTQ files`. In batches 18-19 all three cases
  (`SRR25734143`, `SRR26167346`, `SRR26167354`) were covered by the dataset's duplicate row
  downloading the same file successfully, so nothing was lost. In a non-duplicated dataset the
  run would be dropped and the sample published from partial input. **Batch 22 is that case**,
  and with deduplicated input from batch 22 on, no duplicate row will rescue it again:
  `SRR27373658_2` hit a 16-minute ENA outage (an 802-byte Apache directory listing served for
  the file URL on all 5 attempts, 2026-10-03 11:36-11:52; a valid 3.34 GB gzip the next day),
  `_1` went on alone, and `GSM7996301` was published from 7 of 8 runs, 7.7% of its reads
  missing. GitHub #16.
* **ENA can serve FASTQs that are corrupt at source.** Batch 22, GSE249894: all four files of
  `SRR27178495` and `SRR27178635` fail `pigz -t` (`incomplete deflate data`) on every clean
  download, and each matches ENA's published `fastq_bytes` **and `fastq_md5`** — ENA checksummed
  the truncated files, so MD5 verification (#23) would pass them. `WGET10X` refuses them
  correctly, but the route is fixed as `ENAFQ` in `parse_metadata.sh` and nothing falls back to
  SRA, so the runs are dropped and the samples published partial (`GSM7966587` 6.7%,
  `GSM7966594` 1.8% of reads). **Detection signature**: a `WGET10X` integrity failure on a
  download whose size equals ENA's `fastq_bytes`. Drafted in
  `issues/2026-10-04-batch22-ena-fastq-corrupt-at-source.md`; the fallback it proposes would fix
  #9 too.
* **SRA runs whose read layout changes partway through are refused, with the data complete.**
  `SRA2FASTQ` dumps `--split-files`, which assigns a spot's *n*th read to `_n` regardless of
  role. Batch 22: `SRR27016858`/`SRR27016862` (GSE249159) mix lanes of 3 reads per spot
  (8 + 28 + 91) with lanes of 2 (28 + 91), so `_3` holds about half the spot count and the
  mate-count guard fires; `SRR26923939` (GSE248489) holds R1 and R2 as consecutive single-end
  spots with shared read names and dumps one file. 3 samples lost, all recoverable by splitting
  on read length and pairing on read name. **Detection signature**: `sra.tsv`'s mean spot length
  is not a sum of plausible read lengths (123 = midway between 127 and 119; 58 for 28 + 90 as
  separate spots). Batch 16's `SRR24937677` `partial-mate` may be the same shape. Also: these
  runs retry to 5 on `SRA2FASTQ`'s own `errorStrategy`, 20 deterministic attempts and 59 core-h.
  Drafted in `issues/2026-10-04-batch22-sra-mixed-spot-layout.md`.
* **Nothing gates a sample on whether it mapped.** Batch 9 published 38 near-empty matrices —
  HTO/hashing, ADT/CSP, CRISPR-enrichment and "custom" feature libraries aligned as GEX. Every
  task completed, so none of it is in `failures9.tsv`; the pathology is only visible in
  `mapping_qc_stats.tsv`, and only if you compute it. Neither metadata nor geometry can filter
  these: all 38 are `library_strategy = RNA-Seq` / `transcriptomic single cell`, and GSE216914's
  CSP and VDJ libraries share an identical `10/10/26/90` layout while mapping at 0.9% and 93.9%.
  The wrapper's own strand test already measures the discriminating number and warns without
  acting — a gate at `max(fwd,rev) < 5%` flags 31 samples with zero false alarms, and post hoc
  `exon_u < 0.05` separates the set completely (flagged max 0.0398, kept min 0.0597).
  **Batch 10 validated both halves of that on a second, worse run**: 75 of 748 (10.0%, 967
  CPU-h), and the 5% gate flags 71 of them with zero false alarms — 102 of 113 across 1532
  samples over the two batches, still no false positive. The screen needed one correction, found
  by using it: pair `exon_u < 0.05` with `Med_nFeature < 100`, or single-nucleus data is called
  junk (`GSM6729514`: exon_u 0.043, 583 features, full_u 0.47). `full_u` is **not** the right
  second term — it clears 7 of batch 9's confirmed hashing libraries. The gate belongs in the
  container's `starsolo_10x_auto.sh`; a repo-side `Flag` column for `mapping_qc_stats.tsv` is
  drafted, with diff, in `issues/2026-09-18-batch9-no-mapping-rate-gate.md`. **Do not chase this with a read-length
  threshold** — 42 bp of feature library is junk and 55 bp of GEX maps at 94%, measured in the
  same batch.
* **The STARsolo strand test decides paired-end mode on no margin, and the paired-end STAR call
  segfaults.** `starsolo_10x_auto.sh` sets `PAIRED=True` from `R1LEN > 50` alone, then keeps it
  unless the strand test returns `Forward`; the test concludes `Reverse` on `PCTREV > PCTFWD`,
  strictly greater, with no floor. Samples that stay paired-end run
  `--soloBarcodeMate 1 --clip5pNbases 39 0` and STAR 2.7.10a_alpha_220818 segfaults immediately
  after thread creation. Batch 9's GSE214532: 12 samples, identical 150+150 geometry, 9 lost at
  exit 139 on 6 attempts each (54 attempts, 18 CPU-h), 3 completed at 94-95% mapping — one of
  them on a 44/44 tie. All 12 tripped the script's own "Low percentage of reads mapping to
  GeneFull" warning, so the test was inconclusive for every one of them and was still allowed to
  decide. The file is inside `quay.io/cellgeni/starsolo:v4.3` and there is no flag to override it
  from here; drafted in `issues/2026-09-18-batch9-starsolo-strand-test-segfault.md`. Batch 6's 2
  `star-segfault` rows and batch 3's 1 may be the same branch — unverifiable, their work dirs are
  gone. **Batch 20 cost 52 samples to it** (GSE245310, snRNA-seq, 3 attempts each, ~600 CPU-h),
  and showed that the margin is not the problem: the GeneFull margins were 9-23 points, with
  27 of the 52 at `%rev` ≥ 50. The test reads `GeneFull`, which on single-nucleus data is mostly
  intronic. The test alignments' own exonic `Gene` rows say Forward by 2.6-4.1× in all 52.
  `--wl gex_3pv3_family` had already fixed the chemistry as 3'-only. Skipping the test for a
  3'-only `--wl` saves every case seen. 68 samples across batches 9, 14, 15 and 20.
* **`umi_len` is held to exact equality across a sample's runs, and it is an observed quantity.**
  A run sequenced with a longer R1 reports 12 observed UMI bases where its siblings report 10, and
  `require_metadata_compatibility` rejects the sample. Batch 9 lost GSM6634355 and GSM6634358
  (GSE215362, 9 runs each, 8-vs-1) this way — same chemistry id, same CB length, whitelist match
  93.7% and 95.8%. The run-level inference already tolerates the short case, warning "treating as
  truncated-UMI layout", and then the sample-level guard rejects exactly what it accepted. Fix is
  to take the minimum; drafted in `issues/2026-09-18-batch9-umi-len-mismatch-rejects-sample.md`
  (no diff — registry module with a `.module-info` checksum).
* **Chemistry is still detected twice.** `chemistry.json` reaches `RENAME10XSAMPLE`, but
  `starsolo_10x_auto.sh` re-detects from scratch and can disagree with run-level inference. This
  is `chem-no-whitelist-rand`: 13 late failures in August 2026, still 3 in batch 6.
* **Empty `chemistry.json` reaches `RENAME10XSAMPLE`.** Batch 5's two exit-1
  `RENAME10XSAMPLE` failures are a bare `JSONDecodeError` out of `load_run_metadata` — the file
  exists and is empty. The classifier calls this `chem-json-unreadable`; there is no guard.
* **`FETCH10XMETA` can fail per-sample with no run ID.** `No experiment or run ID found for
  GSM… in GSE….sra.tsv` (batch 6, GSE206528). Not in the August taxonomy at all.
* **No MD5 verification.** ENA publishes `fastq_md5`, but plumbing it through means extending
  `curl_ena_metadata.sh`, the awk join in `parse_metadata.sh`, `collect_metadata.sh` and the
  searchlist schema — and it only covers the ENAFQ route. `gzip -t` catches truncation but not a
  valid gzip holding malformed FASTQ.
* **No pre-download 10x screen.** GSE109816 (880 runs, not 10x) was downloaded in full before
  being rejected at `RENAME10XRUN`. A per-dataset whitelist sample right after `FETCH10XMETA`
  would stop that.
* **Feature-barcode libraries are not filtered on metadata.** ADT/HTO/hashing libraries are
  usually labelled as such in GEO and could be skipped before download.
* **`REPROCESS10X_BAM2FASTQ` stub is broken** — it touches `<id>.fastq.gz` but declares
  `path("fastqs/*")`, blocking stub-testing of any BAM-origin path.
* **`modules/local/reprocess10x/sra2fastq/` is dead code** — the more careful implementation,
  not wired up. Wire it up or delete it.
* **`scripts/run_reprocess.bsub` does not pass `-name`.** Every run gets a generated name and
  the batch → run link has to be recovered from `.nextflow/history`. One flag would fix it.
* **Module drift.** The bsub loads `cellgen/nextflow/26.04.1`, CLAUDE.md and README say
  `26.04.6`, the manifest requires `>=26.04.1`, and the PATH default is 25.04.4 — which reads a
  26.x `.nextflow/cache` as an empty run and rejects the `accelerator` trace field.
