---
title: Outputs
parent: Reference
nav_order: 3
---

# Outputs
{: .no_toc }

1. TOC
{:toc}

## Layout

Everything goes under `--outdir` (default `results/`). Outputs are hard-linked from the work
directory (`--publish_mode link`), so they cost no extra space until the work directory is
cleaned. With `--no_raw`, `raw/` and its index CSVs are not written; everything else is.

```
results/
├── raw/<dataset_id>/
│   ├── fastq/<sample_id>/<run>/      FASTQs as downloaded or dumped
│   ├── bam/<sample_id>/              10x BAM files
│   ├── sra/<sample_id>/              .sra archives
│   ├── runs/<sample_id>/<run>/       per-run FASTQs in Cell Ranger naming
│   └── renamed_fastq/<sample_id>/    the sample's merged FASTQs, as the aligners read them
├── starsolo/<dataset_id>/
│   ├── <sample_id>/                  STARsolo output for one sample
│   └── <dataset_id>.solo_qc.tsv      mapping QC for every sample of the dataset
├── cellranger/<dataset_id>/<sample_id>/   Cell Ranger output (with --cellranger)
├── metadata/<dataset_id>/            what the archives said about the series
├── index/                            CSV indexes of everything above
├── multiqc/                          the run report
├── cleanup/                          with -profile lineage_keep only
└── *.tsv, versions.yml               run-level tables
```

## STARsolo

Each `starsolo/<dataset_id>/<sample_id>/` holds STAR's logs (`Log.final.out`, `Log.out`),
`SJ.out.tab`, `strand.txt` (the strand STARsolo settled on) and `output/`, with one folder per
feature type — `Gene`, `GeneFull` and `Velocyto` — each with `raw/` and `filtered/` matrices,
`Summary.csv` and `Features.stats`. BAMs are not written (`--no-bam`).

## Metadata

`metadata/<dataset_id>/` keeps the raw and parsed archive metadata: `links.tsv` (one row per run:
accession, species, download location, type, sample), the SRA and ENA tables, the sample-to-run
map (`*.sample_x_run.tsv`), the accession lists, and for GEO the SOFT family file or for
ArrayExpress the SDRF/IDF. `links.tsv` is the one to read first.

## Indexes

Each CSV in `index/` lists one kind of published output, one row per item, with its identifiers
(run, sample, dataset as applicable) and the path:

| File | Lists |
|---|---|
| `starsolo.csv` | STARsolo outputs, with the species and the chemistry passed to each aligner |
| `cellranger.csv` | Cell Ranger outputs, likewise |
| `qc.csv` | The per-dataset `solo_qc.tsv` files |
| `original_fastq.csv` | FASTQs as downloaded or dumped |
| `runs.csv` | Per-run renamed FASTQs |
| `renamed_fastq.csv` | Per-sample merged FASTQs |
| `bam.csv`, `sra.csv` | BAMs and `.sra` archives |

The aligner indexes always have the same columns (`id`, `dataset_id`, `specie`, `chemistry`,
`wl`, `path`); `chemistry` and `wl` are empty for a sample whose chemistry was left to the
aligner.

## Run-level tables

Written by REPROCESS10X_COLLATE, one of the last two tasks of a run:

| File | One row per | Columns |
|---|---|---|
| `datasets.tsv` | dataset | samples requested (all = the whole series), found in the metadata, resolved, skipped and mapped; runs; data types; species |
| `samples.tsv` | sample | resolved species, Cell Ranger chemistry, STARsolo whitelist id, number of runs, origin (the runs' download types — `ENAFQ`, `ORIFQ`, `SRA`, `BAM`, or several), whether it was mapped |
| `runs.tsv` | run | sample, dataset, download type, species |
| `run_chemistry.tsv` | run | inferred chemistry, layout confidence, modality, whitelist, barcode and UMI lengths, read lengths, strand hint, warnings |
| `mapping_qc_stats.tsv` | mapped sample | reads, reads and fraction in cells, UMIs in cells, cells, median features, good-barcode fraction, whitelist, species, paired, strand, unique and multi-mapping rates (all, exon, full gene) |
| `skipped_samples.tsv` | skipped sample | samples dropped for want of a usable species (absent if none) |
| `versions.yml` | process | software versions |

## Run report and cleanup

`multiqc/multiqc_report.html` and `multiqc/multiqc_data/` — see [Run report](run-report.md).

With `-profile lineage_keep`, `cleanup/` holds `deleted.tsv`, `processes.tsv`,
`failed_tasks.tsv` and `rerun.tsv` — see [Reruns and cleanup](../guides/reruns.md).

Nextflow's own timeline, execution report, trace and DAG go to `reports/` in the launch
directory, not under `--outdir`.
