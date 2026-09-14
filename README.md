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
| `nf-test.config` | nf-test configuration |
| `tests/config/nf-test.config` | Nextflow config applied to test runs |
| `tests/scripts/run_tests.bsub` | LSF submission script for the test suite |

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
| `--starsolo` | Run STARsolo alignment after downloading | `false` |
| `--cellranger` | Run Cell Ranger alignment after downloading (not yet implemented) | `false` |
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
├── raw/
│   └── <dataset_id>/
│       ├── fastq/
│       │   └── <sample_id>/        FASTQs (R1, R2, I1)
│       ├── bam/
│       │   └── <sample_id>/        10x BAM files
│       └── sra/
│           └── <sample_id>/        SRA files
├── starsolo/
│   └── <dataset_id>/               STARsolo count matrices and QC stats
├── metadata/
│   └── <dataset_id>/               Metadata files (links, parsed TSVs, SOFT, etc.)
├── index/
│   ├── fastq.csv                   Index of all published FASTQs
│   ├── bam.csv                     Index of all published BAMs
│   ├── sra.csv                     Index of all published SRA files
│   └── starsolo.csv                Index of all STARsolo outputs
├── versions.yml                    Software versions used by each process
└── mapping_qc_stats.tsv            Per-sample STARsolo mapping QC statistics
```

## Species handling

By default (`--no_infer_specie` not set), species is read from sample metadata. If it is blank or unrecognised:

- If `--default_specie` is set, that species is used and a warning is logged.
- If `--default_specie` is not set, the sample is assigned `UNKNOWN` and skipped by STARsolo.

Use `--no_infer_specie` to bypass metadata entirely and force all samples to `--default_specie`.

## Tests

The test suite uses [nf-test](https://www.nf-test.com). Submit it to LSF:

```bash
bsub < tests/scripts/run_tests.bsub
```

The job is only the driver — Nextflow submits the test processes themselves to the
`transfer` queue. Results land in `logs/nf-test/testOutput.<jobid>.log`.

Pass extra nf-test arguments through `ARGS`:

```bash
# just the regression tests for GEO series with missing sample relations
bsub -env "all, ARGS=--tag relation-recovery" < tests/scripts/run_tests.bsub

# re-record the snapshots after an intentional change
bsub -env "all, ARGS=--update-snapshot" < tests/scripts/run_tests.bsub

# fail on a snapshot mismatch instead of silently updating it
bsub -env "all, ARGS=--ci" < tests/scripts/run_tests.bsub
```

Available tags: `geo`, `arrayexpress`, `bioproject`, `enafq`, `orifq`, `bam`, `sra`,
`mouse`, `subset`, `relation-recovery`, `stub`.

The full suite takes around ten minutes. See
[modules/cellgeni/fetch10xmeta/README.md](modules/cellgeni/fetch10xmeta/README.md#tests)
for what each test covers.

### Running outside LSF

The tests download live metadata from GEO, SRA, ENA and BioStudies, so they need
outbound network access — a farm head node or the `transfer` queue. To run them
directly:

```bash
module load cellgen/nf-test/0.9.5
module load cellgen/nextflow/26.04.6     # the manifest requires >=26.04.1

cd modules
nf-test test cellgeni/fetch10xmeta/tests/main.nf.test --config ../nf-test.config --profile local
```

Two things about that invocation are deliberate:

- **`cd modules`.** Before running anything, nf-test walks the entire launch
  directory to build a dependency graph. It does this regardless of `testsDir`,
  and the `ignore` option only filters the result afterwards. Launched from the
  repository root that walk never returns, because `nf-work` holds hundreds of
  thousands of files on Lustre. `modules` holds under a hundred and is the only
  tree the tests need; snapshots are still written next to each test.
- **`--profile local`** runs the test processes on the current host, where
  `--profile lsf` submits them to the `transfer` queue.

## Requirements

- Nextflow `>=26.04.1`
- Singularity (or Docker for local runs)
- LSF cluster (or adjust `nextflow.config` executor for local use)
- For the tests: nf-test `>=0.9.2` and the [nft-csv](https://github.com/lukfor/nft-csv)
  plugin, which nf-test downloads on first run
