---
title: Datasets table
parent: Reference
nav_order: 2
---

# Datasets table
{: .no_toc }

1. TOC
{:toc}

## Format

`--datasets` takes a **tab-separated** file, whatever its extension, with a header row and two
columns:

```
dataset_id	sample_id
GSE230685	GSM7232572,GSM7232573
E-MTAB-9221	ERS4689152,ERS4689153
PRJEB37166	ERS4605100,ERS4605101
```

| Column | Holds |
|---|---|
| `dataset_id` | A GEO series (`GSE*`), BioProject (`PRJNA*`, `PRJEB*`, `PRJDB*`) or ArrayExpress accession (`E-MTAB-*`) |
| `sample_id` | Comma-separated sample accessions in that dataset: `GSM*` for GEO, otherwise the SRA sample accessions (`SRS*`, `ERS*`, `DRS*`). Empty for the whole series |

The commas inside `sample_id` are why the file must be tab-separated: a `.csv` extension does not
change how it is read. See [`examples/datasets.tsv`](https://github.com/cellgeni/nf-reprocessing-public-10x/blob/main/examples/datasets.tsv)
for a longer example.

## An empty `sample_id` means the whole series

{: .warning }
Leaving `sample_id` empty processes **every** sample in the dataset. That is how you ask for a
whole series, but it is also easy to do by accident: a table that declares the header and
supplies only `dataset_id` values looks correct on inspection. A series is often mixed, so the
cost can be large — one ticket fetched ~596.8 GB of bulk RNA-seq, scATAC and bisulfite data this
way before any of it was rejected downstream.

Run with `--metaonly` first to see what a series actually contains before committing the
download: `results/metadata/<dataset_id>/links.tsv` lists every run found, with its species, data
type and location.

## What the archives return

For each dataset, FETCH10XMETA reads the archive metadata (NCBI SRA, EBI ENA or BioStudies,
depending on the accession), resolves each sample to its runs, and picks where to download each
run from:

| Type | Downloaded as | Then |
|---|---|---|
| `ORIFQ` | The submitter's original FASTQs | renamed |
| `ENAFQ` | ENA's FASTQs | renamed |
| `BAM` | A 10x BAM | converted back to FASTQ (`bamtofastq`) |
| `SRA` | An `.sra` archive | checked, then dumped to FASTQ; archives the dump tools cannot read (a NULL `NAME` column or a `GenericFastq` schema) are rebuilt with `vdb-dump` instead |

A sample can mix types across its runs; they are merged before renaming, so each sample gives one
set of FASTQs and one aligner job. `results/runs.tsv` and `results/datasets.tsv` record what was
found for each run and dataset.
