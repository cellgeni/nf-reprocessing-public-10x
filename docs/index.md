---
title: Home
nav_order: 1
---

# nf-reprocessing-public-10x

nf-reprocessing-public-10x fetches public 10x single-cell datasets from GEO, SRA, ENA or
ArrayExpress and reprocesses them uniformly. Give it a table of series and samples; it reads the
archive metadata, downloads whatever the submitters deposited (FASTQ, 10x BAM or `.sra`),
normalises it to Cell Ranger FASTQ naming, infers the chemistry from the reads themselves, and
aligns with STARsolo or Cell Ranger. One run report covers every sample.

```
datasets.tsv ──► archive metadata ──► FASTQ / BAM / SRA ──► Cell Ranger naming ──► chemistry
                 (FETCH10XMETA)        (download, convert)    (per run, per sample)  inference
                                                                                        │
                       run report ◄── run tables ◄── STARsolo / Cell Ranger ◄──────────┘
                       (MultiQC)      (COLLATE)       (per species)
```

## Quick start

On the Sanger farm:

```bash
module load cellgen/nextflow/26.04.6

nextflow run main.nf -profile sanger \
  --datasets datasets.tsv --default_specie human --starsolo
```

Anywhere else, choose a container engine and pass the references, which only `-profile sanger`
provides (see [Profiles and running elsewhere](guides/profiles.md)):

```bash
nextflow run main.nf -profile docker \
  --datasets datasets.tsv --default_specie human --starsolo \
  --wl_basedir whitelists --star_human_reference /path/to/star/GRCh38
```

`datasets.tsv` is tab-separated, one row per series:

```
dataset_id	sample_id
GSE230685	GSM7232572,GSM7232573
PRJEB37166	ERS4605100,ERS4605101
```

Results go to `results/`: STARsolo output per dataset in `results/starsolo/`, run tables at the
top level and the run report in `results/multiqc/`.

{: .warning }
An empty `sample_id` processes the **whole series**. Run with `--metaonly` first to see what a
series holds — see [Datasets table](reference/datasets.md).

## Where next

| Page | Read it when |
|---|---|
| [Getting started](getting-started.md) | You are running the pipeline for the first time |
| [Profiles and running elsewhere](guides/profiles.md) | You are choosing a profile, or running outside the Sanger farm |
| [Resources and queues](guides/resources.md) | You want to know what each task asks for, or tune it for your machine |
| [Species](guides/species.md) | Samples were skipped for want of a species, or the metadata names the wrong one |
| [Chemistry inference](guides/chemistry.md) | You want to know how the chemistry is decided, or why a run was rejected |
| [Reruns and cleanup](guides/reruns.md) | A run stopped part-way, or you want the work dirs cleaned as it goes |
| [Parameters](reference/parameters.md) | You need the full list of parameters and what each accepts |
| [Outputs](reference/outputs.md) | You want to know what each output file is |
| [Run report](reference/run-report.md) | You are reading the MultiQC report and want to know where to look |
| [Troubleshooting](troubleshooting.md) | Something failed with an error message |
| [Knowledge base](knowledge-base.md) | A public dataset misbehaved and you want to know whether it is a known problem |
| [Developers](developers.md) | You are changing the pipeline or its modules |
| [Changelog](changelog.md) | You want to know what changed between releases |
