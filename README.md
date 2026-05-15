# nf-reprocessing-public-10x

Nextflow pipeline for loading and reprocessing public 10x datasets from GEO, SRA, ENA, or ArrayExpress.

## Repo structure

| Path | Description |
|---|---|
| `main.nf` | Main Nextflow pipeline entry point |
| `nextflow.config` | Pipeline configuration — LSF executor, Singularity, and default params |
| `workflow/main.nf` | Core workflow: metadata fetch, download, and STARsolo alignment |
| `subworkflows/` | Download and STARsolo alignment subworkflows |
| `modules/` | Individual process modules |
| `examples/datasets.tsv` | Example input file |
| `examples/RESUME` | Example run script |

## Usage

```bash
nextflow run main.nf --datasets <datasets.tsv> [OPTIONS]
```

### Parameters

| Parameter | Description | Default |
|---|---|---|
| `--datasets` | Path to a TSV file with dataset and sample IDs | required |
| `--outdir` | Directory to save results | `results` |
| `--default_specie` | Species to assign when metadata is missing or unknown (`human` or `mouse`). Samples without a resolved species are skipped by STARsolo. | `null` |
| `--no_infer_specie` | Skip reading species from metadata; assign `--default_specie` to all samples. Requires `--default_specie`. | `false` |
| `--metaonly` | Only fetch metadata — skip downloading and alignment | `false` |
| `--help` | Print help message and exit | — |

### Input file format

The `--datasets` TSV must have a header row with two columns:

```tsv
dataset_id	sample_id
GSE230685	GSM7232572,GSM7232573
E-MTAB-9221	ERS4689152,ERS4689153
PRJEB37166	ERS4605100,ERS4605101
```

- `dataset_id` — GEO series (`GSE*`), BioProject (`PRJEB*`/`PRJNA*`), or ArrayExpress accession (`E-MTAB-*`)
- `sample_id` — comma-separated list of sample accessions belonging to that dataset

See [examples/datasets.tsv](examples/datasets.tsv) for a full example.

## Quick example

```bash
nextflow run main.nf \
  --datasets examples/datasets.tsv \
  --outdir results \
  --default_specie human \
  -resume
```

## Output structure

```
results/
  raw/<dataset_id>/fastq/<sample_id>/   FASTQs
  raw/<dataset_id>/bam/<sample_id>/     10x BAM files
  raw/<dataset_id>/sra/<sample_id>/     SRA files
  starsolo/<dataset_id>/                STARsolo count matrices and QC stats
  metadata/<dataset_id>/                Metadata files (links, parsed TSVs, SOFT, etc.)
  index/fastq.csv                       Index of all published FASTQs
  index/bam.csv                         Index of all published BAMs
  index/sra.csv                         Index of all published SRA files
  index/starsolo.csv                    Index of all STARsolo outputs
  versions.yml                          Software versions used by each process
  mapping_qc_stats.tsv                  Per-sample STARsolo mapping QC statistics
```

## Species handling

By default (`--no_infer_specie` not set), species is read from sample metadata. If it is blank or unrecognised:

- If `--default_specie` is set, that species is used and a warning is logged.
- If `--default_specie` is not set, the sample is assigned `UNKNOWN` and skipped by STARsolo.

Use `--no_infer_specie` to bypass metadata entirely and force all samples to `--default_specie`.

## Requirements

- Nextflow `>=26.04.1`
- Singularity (or Docker for local runs)
- LSF cluster (or adjust `nextflow.config` executor for local use)
