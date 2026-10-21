---
title: Changelog
nav_order: 8
---

# Changelog
{: .no_toc }

All notable changes to nf-reprocessing-public-10x. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Releases are tagged `YY-DDD`: the year
and the day of the year. The [GitHub releases](https://github.com/cellgeni/nf-reprocessing-public-10x/releases)
have the full notes for each.

1. TOC
{:toc}

## Unreleased

### Added
- **`docker`, `singularity` and `sanger` profiles.** The pipeline runs outside the Sanger farm:
  the base configuration names no executor, queue or container engine, and `-profile sanger`
  carries everything farm-specific — LSF, the Singularity cache and binds, the farm's whitelists
  and references, and the queues. See [Profiles and running elsewhere](guides/profiles.md).
- **Resource labels**, as in nf-core: `process_single`, `process_high`, `process_long`, and
  `process_network` for the tasks that fetch from the archives. `configs/sanger.config` picks the
  queue from each task's resolved memory and time instead of naming one per process.
- **Launch checks** in `main.nf`: a run with no container engine, without `--wl_basedir`
  (unless `--metaonly`), or with `--starsolo`/`--cellranger` and no reference for them now stops
  at once. Each used to exit 0 having done nothing.
- **This documentation site**, on Just the Docs: getting started, guides, reference (with the
  parameters generated from `nextflow_schema.json`), troubleshooting, developers and this
  changelog, beside the existing knowledge base. The README is now an overview.

### Changed
- **`-profile sanger` is required on the farm**, on every run and `-resume`
  (`-profile sanger,lineage_keep` for the batch scripts). Without a profile the pipeline refuses to
  start.
- Every task now asks for an explicit run limit. The values match the limits of the queues the
  tasks used before (12 h, 48 h), so the farm routes exactly as it did.
- The reference parameters have no defaults outside `-profile sanger`.

## 26-292 — 2026-10-19

Resource prediction, SRA fixing, run reports, `lineage_keep`.

### Added
- **STARsolo resources planned per sample.** `lib/StarsoloResources.groovy` estimates reads from
  the head of each R2 FASTQ and sets memory, run limit and queue for each attempt; the
  coefficients are fitted in `notebooks/resource_prediction.ipynb`.
- **Every `.sra` is checked before it is dumped** (`REPROCESS10X_SRACHECK`). Runs with a NULL
  `NAME` column or a `GenericFastq` schema, which `fastq-dump` and `fasterq-dump` cannot read, are
  rebuilt with `vdb-dump` by the new `REPROCESS10X_SRAFIX`.
- **Run-level tables and a MultiQC run report** (`REPROCESS10X_COLLATE`, `REPROCESS10X_RUNREPORT`),
  including *Samples without a mapped output*, the run's failure manifest.
- **`-profile lineage_keep`**: the `plugins/nf-lineage-keep` plugin empties work directories as
  samples finish and keeps everything a failure needs; `cleanup/rerun.tsv` lists the failed
  samples for a `-resume`.
- **nf-schema validation** (`nextflow_schema.json`), so `--help` lists every parameter.
- **`--no_raw`**, which stops publishing `raw/` so a cleaned run actually frees its FASTQs.
- Knowledge base: post-mortems for batches 7-25, archive pathologies through batch 24, and the
  [tchem benchmark](tchem-benchmark.md).

### Changed
- `SRA2FASTQ` asks for 16 GB first and 128 GB on retry.
- STARsolo container `v4.4`.
- Nextflow 26.04.6 or later is required.

### Fixed
- An `SRACHECK` that was ignored after its retries crashed the join that feeds `SRA2FASTQ`.
- The failure reporters omitted `SRAFIX`; the run report invented a sample for each BAM
  conversion failure.

{: .note }
Adding `meta.resources` and moving to STARsolo v4.4 change STARsolo's task hash, so a batch
started before this release re-runs STARsolo on `-resume`.

## 26-257 — 2026-09-14

STARsolo or Cell Ranger, inferred chemistry, guarded downloads.

### Added
- **Two aligners, one path.** `--starsolo` and `--cellranger` select the aligner; both can run in
  one execution, each per species.
- **Chemistry inferred from the reads**, once per sample, and translated separately for STARsolo's
  `--wl` and Cell Ranger's `--chemistry`. Inference samples a window of reads rather than the head
  of the file. The [10x chemistry reference](10x_chemistry_reference.md) documents the table it
  implements.
- **`--metaonly`**, to inspect a series before committing to the download.
- **`skipped_samples.tsv`** for samples dropped for want of a usable species.
- The documentation knowledge base and the `reprocess-debug` skill; an nf-test suite for
  FETCH10XMETA.

### Changed
- **`links.tsv` is read whole**, so species, run counts and alignable-sample counts are settled
  before they become grouping keys.
- BAM-derived and FASTQ/SRA-derived runs of one sample merge before renaming: one renamed set and
  one aligner job per sample. Mapping QC runs once per dataset.
- Metadata parsing consolidated into one `parse_metadata.sh`, with sample-to-run relations
  recovered when the subset file is missing.

### Fixed
- `WGET10X` rejects empty URLs, expects exactly one file per task, and tests the archive.
- `SRA2FASTQ` refuses to emit a partial run: single-mate dumps, mates that disagree and zero-length
  reads fail at the guard instead of two stages downstream.
- I1/I2 index reads are no longer dropped between `RENAME10XRUN` and `RENAME10XSAMPLE`.
- Boolean flags are validated: `--starsolo ture` used to parse as false, download a whole batch,
  align nothing and exit 0.

## 26-138 — 2026-05-18

- The metadata fetch moved to the `cellgeni/fetch10xmeta` registry module.
- `--metaonly` and `--no_infer_specie` validated as booleans.

## 26-135 — 2026-05-15

- The pipeline restructured: a separate workflow for the pipeline logic, input validation, and
  publishing through the workflow `output {}` block instead of per-process `publishDir`.
- Metadata parsing and SRA conversion hardened.

## 23-263 — 2023-09-20

The first tagged release: three-step reprocessing (metadata, download and conversion to Cell
Ranger FASTQ naming, STARsolo), with matching on sample ids throughout.
