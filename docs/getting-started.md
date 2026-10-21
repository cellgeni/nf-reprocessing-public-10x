---
title: Getting started
nav_order: 2
---

# Getting started
{: .no_toc }

1. TOC
{:toc}

## Requirements

- **Nextflow 26.04.6 or later**; the manifest refuses anything older. On the Sanger farm the
  default on `PATH` is older, so `module load cellgen/nextflow/26.04.6` first.
- **A container engine**: Singularity or Docker. Every step except the metadata fetch runs in a
  public image from `quay.io`; the pipeline stops at launch if no engine is enabled.
- **Outbound network access** for the tasks that fetch metadata and data.
- **Disk.** Raw data is large: a single series can run to hundreds of GB of FASTQ. Check
  `--metaonly` output before committing (step 2).
- **On the Sanger farm**: LSF, and membership of the `cellgeni` group for the shared image cache
  and references. Elsewhere: a STAR index and the barcode whitelists — see
  [Profiles and running elsewhere](guides/profiles.md#running-outside-the-sanger-farm).

The [nf-schema](https://registry.nextflow.io/plugins/nf-schema@3.0.0) plugin, which checks the
parameters, is downloaded by Nextflow on the first run.

## 1. Write a datasets table

A tab-separated file, whatever its extension, with a header and one row per series:

```
dataset_id	sample_id
GSE230685	GSM7232572,GSM7232573
E-MTAB-9221	ERS4689152,ERS4689153
PRJEB37166	ERS4605100,ERS4605101
```

`dataset_id` is a GEO series, BioProject or ArrayExpress accession; `sample_id` is a
comma-separated list of its samples. [Datasets table](reference/datasets.md) has the details.

{: .warning }
Leaving `sample_id` empty processes the entire series. That is how you ask for a whole dataset,
but series are often mixed: one ticket fetched ~597 GB of bulk RNA-seq, scATAC and bisulfite
data this way before any of it was rejected.

## 2. Look before you download

```bash
nextflow run main.nf -profile sanger --datasets datasets.tsv --metaonly
```

`--metaonly` fetches only the archive metadata. `results/metadata/<dataset_id>/links.tsv` lists
every run it found, with its species, data type and download location, and
`results/datasets.tsv` counts samples requested, found and resolved per series.

## 3. Run it

### On the Sanger farm

The Nextflow driver must stay alive until every task has finished, so submit it as an LSF job of
its own rather than running it on a head node:

```bash
#BSUB -G cellgeni
#BSUB -q oversubscribed
#BSUB -n 1
#BSUB -M 4GB
#BSUB -R "select[mem>4GB] rusage[mem=4GB]"
#BSUB -o logs/reprocess.%J.log
#BSUB -e logs/reprocess.%J.err

module load cellgen/nextflow/26.04.6

nextflow run main.nf -profile sanger \
    --datasets datasets.tsv \
    --outdir results \
    --default_specie human \
    --starsolo \
    -ansi-log false
```

`-profile sanger` sets LSF, the group's Singularity cache, the farm's whitelists and references,
and the queues. Add `lineage_keep` (`-profile sanger,lineage_keep`) to clean work dirs as samples
finish — see [Reruns and cleanup](guides/reruns.md).

### Anywhere else

```bash
nextflow run main.nf -profile docker \
    --datasets datasets.tsv --default_specie human --starsolo \
    --wl_basedir whitelists \
    --star_human_reference /path/to/star/GRCh38 \
    -c local.config
```

[Profiles and running elsewhere](guides/profiles.md) explains where the whitelists come from,
what goes in `local.config`, and how to point the pipeline at a cluster.

## 4. Check the results

```
results/
├── starsolo/<dataset_id>/     STARsolo output per sample, and the dataset's QC table
├── metadata/<dataset_id>/     what the archives said about each series
├── raw/<dataset_id>/          downloads and FASTQs (not with --no_raw)
├── index/                     CSV indexes of everything published
├── multiqc/                   the run report
└── datasets.tsv, samples.tsv, runs.tsv, run_chemistry.tsv, mapping_qc_stats.tsv, versions.yml
```

Start with `results/multiqc/multiqc_report.html`. Its **Samples without a mapped output** table
lists every sample that did not reach a STARsolo output, with the stage and exit class it
stopped at. [Run report](reference/run-report.md) explains each section and
[Outputs](reference/outputs.md) each file.

{: .important }
A task that fails twice is ignored, not fatal, so a run can finish with exit status 0 and still
have lost samples. The run report is where they show up.

## Getting help on the command line

```bash
nextflow run main.nf --help
```

prints every parameter with its type and default, generated from `nextflow_schema.json` — the
same list as [Parameters](reference/parameters.md).

## When a run is refused

Parameters are checked before anything is fetched. A wrong type, an unknown species or a missing
`--datasets` stops the run with nf-schema's list of invalid values; the value in brackets is
what it got. After that, `main.nf` checks what the schema cannot: a container engine, and the
whitelists and references the chosen aligners need. Each of these stops the run with a message
saying what to add — see [Troubleshooting](troubleshooting.md#the-run-stops-at-launch).
