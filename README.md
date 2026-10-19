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
| `lib/StarsoloResources.groovy` | STARsolo memory, run limit and queue predicted from each sample's read count |
| `nextflow_schema.json` | Parameter schema: `--help`, validation and the launch summary (nf-schema) |
| `notebooks/resource_prediction.ipynb` | Rerun/waste analysis and the fit behind `lib/StarsoloResources.groovy` |
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

Parameters are defined in [`nextflow_schema.json`](nextflow_schema.json) and validated at launch
by [nf-schema](https://github.com/nextflow-io/nf-schema): a wrong type (`--starsolo ture`), an
unknown species or a missing `--datasets` stops the run before anything is fetched.
`nextflow run main.nf --help` prints the full list; the main ones:

| Parameter | Description | Default |
|---|---|---|
| `--datasets` | Path to a TSV file with dataset and sample IDs | required |
| `--outdir` | Directory to save results | `results` |
| `--default_specie` | Species to assign when metadata is missing or unknown (`human` or `mouse`). Samples without a resolved species are skipped by STARsolo. | `null` |
| `--no_infer_specie` | Skip reading species from metadata; assign `--default_specie` to all samples. Requires `--default_specie`. | `false` |
| `--metaonly` | Only fetch metadata — skip downloading and alignment | `false` |
| `--starsolo` | Run STARsolo alignment after downloading | `false` |
| `--cellranger` | Run Cell Ranger alignment after downloading (not yet implemented) | `false` |
| `--no_raw` | Publish nothing under `raw/` (downloads, per-run and renamed FASTQs) and none of their index CSVs; the files then live only in the work dir | `false` |
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

**Leaving `sample_id` empty processes the entire series.** This is intentional and is how you
ask for a whole dataset, but it is also easy to do by accident — a table that declares the
header and supplies only `dataset_id` values looks correct on inspection. A series is often
mixed, so the cost can be large: one ticket fetched ~596.8 GB of bulk RNA-seq, scATAC and
bisulfite data this way before any of it was rejected downstream. Run `--metaonly` first to see
what a series actually contains before committing the download.

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

With `--no_raw`, `raw/` and its index CSVs are not written; everything else is.

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
├── multiqc/
│   ├── multiqc_report.html         Run report (see below)
│   └── multiqc_data/               Everything the report plots, as tables
├── mapping_qc_stats.tsv            Per-sample STARsolo mapping QC statistics
├── datasets.tsv                    Per dataset: samples requested / found / resolved / skipped / mapped, runs
├── samples.tsv                     Per sample: resolved species, chemistry, whitelist, runs and origin
├── runs.tsv                        Per run: sample, dataset, data type, species
├── run_chemistry.tsv               Per run: inferred chemistry, confidence, read lengths, warnings
├── skipped_samples.tsv             Samples dropped for want of a usable species (absent if none)
└── versions.yml                    Software versions used by each process
```

The run-level tables are written by `REPROCESS10X_COLLATE` and the report by
`REPROCESS10X_RUNREPORT`, the last two tasks of a run. The report (MultiQC) covers task outcomes
and wasted compute per process, the status of every sample (mapped, skipped, or the stage and exit
code it failed at), chemistry calls, STARsolo mapping metrics, predicted vs observed STARsolo
resources, and software versions.

## STARsolo resources

STARsolo's memory request, run limit and LSF queue are set per sample from its input read count,
which is estimated before the task starts by reading the head of each R2 FASTQ
([`lib/StarsoloResources.groovy`](lib/StarsoloResources.groovy)). Samples whose limit fits in 12 h
go to `normal` (which backfills), the rest to `long`; each retry raises memory ×1.5 and time ×2.
The coefficients come from [`notebooks/resource_prediction.ipynb`](notebooks/resource_prediction.ipynb),
fitted on batches 9-21; refit there and paste the printed constants into the class.

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
