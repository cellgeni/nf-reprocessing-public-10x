# nf-reprocessing-public-10x

Nextflow pipeline for loading and reprocessing public 10x datasets from GEO, SRA, ENA, or ArrayExpress.

## Repo structure

| Path | Description |
|---|---|
| `main.nf` | Main Nextflow pipeline |
| `nextflow.config` | Pipeline configuration — LSF executor, Singularity, and default params |
| `examples/datasets.tsv` | Example input file with dataset and sample IDs |
| `examples/RESUME` | Example run script |
| `subworkflows/` | Download and STARsolo alignment subworkflows |
| `modules/` | Individual process modules |

## Usage

```bash
nextflow run main.nf --datasets <datasets.tsv> [OPTIONS]
```

### Parameters

| Parameter | Description | Default |
|---|---|---|
| `--datasets` | Path to a TSV file with dataset and sample IDs | required |
| `--output_dir` | Directory to save results | `results` |
| `--help` | Print help message and exit | — |

### Input file format

The `--datasets` TSV file must have a header row with two columns:

```tsv
dataset_id	sample_id
GSE230685	GSM7232572,GSM7232573
E-MTAB-9221	ERS4689152,ERS4689153
PRJEB37166	ERS4605100,ERS4605101
```

- `dataset_id` — GEO series (GSE), ENA study (PRJEB/PRJNA), or ArrayExpress accession (E-MTAB-*)
- `sample_id` — comma-separated list of sample accessions belonging to that dataset

See [examples/datasets.tsv](examples/datasets.tsv) for a full example.

## Quick example

```bash
nextflow run main.nf \
  -profile standard \
  --datasets examples/datasets.tsv \
  --output_dir results \
  --ansi-log false \
  -resume
```

## Output

Results are written to `--output_dir` (default: `results/`):

| File | Description |
|---|---|
| `mapping_qc_stats.tsv` | Per-sample STARsolo mapping QC statistics |
| `versions.yml` | Software versions for all pipeline steps |
| `<dataset_id>/` | STARsolo count matrices per dataset |

## Profiles

| Profile | Description |
|---|---|
| `standard` | LSF cluster at Sanger (default) — uses Singularity |
| `development` | Local execution — uses Docker |
