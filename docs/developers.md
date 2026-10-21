---
title: Developers
nav_order: 7
---

# Developers
{: .no_toc }

1. TOC
{:toc}

## Repository layout

| Path | Contents |
|---|---|
| `main.nf` | Entry point: help text, parameter validation and the launch checks, the run tables and run report, and the `output {}` block that decides what is published where |
| `workflow/main.nf` | `REPROCESS10X`, the real orchestrator: metadata, download, the aligners per species, mapping QC |
| `subworkflows/local/download10x` | Download and conversion to per-sample Cell Ranger FASTQs; chemistry resolution |
| `subworkflows/local/starsolo10x` | STARsolo per species, with the per-sample resource plan |
| `modules/cellgeni/` | Modules in Nextflow registry layout (see [Modules](#modules)) |
| `modules/local/reprocess10x/` | Repo-local modules |
| `lib/StarsoloResources.groovy` | STARsolo memory and run limit predicted from each sample's read count |
| `nextflow.config` | Default parameters, error strategy, reports, and the profiles |
| `configs/base.config` | Resources by label |
| `configs/sanger.config` | The `sanger` profile: LSF, Singularity cache, farm references, queue rule |
| `configs/<process>.config` | Per-process exceptions to the labels: `ext.args`, retries, memory rules |
| `nextflow_schema.json` | Parameter types, defaults and descriptions; drives validation, `--help` and [Parameters](reference/parameters.md) |
| `plugins/nf-lineage-keep` | The work-dir cleanup plugin behind `-profile lineage_keep` |
| `notebooks/resource_prediction.ipynb` | Rerun and waste analysis, and the fit behind `StarsoloResources` |
| `docker/` | A Dockerfile for the download and conversion tools (SRA Toolkit, samtools, bamtofastq, STAR, seqtk) |
| `tests/`, `nf-test.config` | nf-test configuration and the LSF script that runs the suite |
| `docs/` | This site |
| `.claude/skills/reprocess-debug/` | The agent skill for debugging finished runs (see [For maintainers](knowledge-base.md#for-maintainers)) |

`main` holds the released pipeline and `dev` the next release; feature branches merge into
`dev`. Releases are tagged `YY-DDD`, the year and the day of the year.

## Architecture

```
main.nf                     nf-schema validation, launch checks, publish:/output: blocks, the two end-of-run tasks
├── workflow/main.nf        REPROCESS10X
│   ├── modules/cellgeni/fetch10xmeta          → links.tsv, one row per run
│   ├── subworkflows/local/download10x         → WGET10X → {BAM2FASTQ, SRACHECK → SRA2FASTQ | SRAFIX} → RENAME10XRUN → RENAME10XSAMPLE
│   ├── subworkflows/local/starsolo10x         → instantiated twice, once per species; sets meta.resources
│   ├── modules/cellgeni/cellranger/count      → likewise, human + mouse
│   └── modules/local/reprocess10x/mappingqc   → one QC job per dataset, across all its species
├── modules/local/reprocess10x/collate         → run-level tables
└── modules/local/reprocess10x/runreport       → parses the run's own trace + those tables, then MultiQC
```

Four ideas hold the wiring together; read the comments where they live before changing them.

- **`links.tsv` is read whole, not row by row.** Species, per-sample run counts and the
  dataset's alignable-sample count are all settled before any becomes a grouping key, because
  none can be decided from one run's row.
- **`meta.sample_id` and `meta.dataset_id` are `GroupKey`s.** They carry the group size the
  downstream `groupTuple` needs; `.getGroupTarget()` reads the accession back, and `main.nf`'s
  `unwrapGroupKeys()` does it before publishing.
- **Chemistry is resolved once per sample and pushed into `meta`**, then translated separately
  for each aligner ([Chemistry inference](guides/chemistry.md)).
- **Runs of mixed origin are merged before RENAME10XSAMPLE**, on a key that deliberately omits
  species, so each sample makes exactly one aligner job.

## Modules

`modules/cellgeni/*` are published to the [Nextflow module registry](https://registry.nextflow.io/);
`.module-info` records the published checksum, and Nextflow warns at launch when a local copy
differs. Two are imported bare (`from 'cellgeni/fetch10xmeta'`, `'cellgeni/starsolo10x'`), the
rest by relative path. `modules/local/*` are repo-local.

| Module | Process | Container |
|---|---|---|
| `cellgeni/fetch10xmeta` | FETCH10XMETA | none: the host's `curl` |
| `cellgeni/wget10x` | WGET10X | `quay.io/cellgeni/reprocess_10x:latest` |
| `cellgeni/sra2fastq` | SRA2FASTQ | `quay.io/cellgeni/reprocess_10x:latest` |
| `local/reprocess10x/sracheck` | REPROCESS10X_SRACHECK | `quay.io/cellgeni/reprocess_10x:latest` |
| `local/reprocess10x/srafix` | REPROCESS10X_SRAFIX | `quay.io/cellgeni/reprocess_10x:latest` |
| `local/reprocess10x/bam2fastq` | REPROCESS10X_BAM2FASTQ | `quay.io/cellgeni/reprocess_10x:latest` |
| `cellgeni/rename10xrun` | RENAME10XRUN | `quay.io/cellgeni/metacells-python:latest` |
| `cellgeni/rename10xsample` | RENAME10XSAMPLE | `quay.io/cellgeni/metacells-python:latest` |
| `cellgeni/starsolo10x` | STARSOLO10X | `quay.io/cellgeni/starsolo:v4.4` |
| `cellgeni/cellranger/count` | CELLRANGER_COUNT | `quay.io/nf-core/cellranger:9.0.1` |
| `local/reprocess10x/mappingqc` | REPROCESS10X_MAPPINGQC | `quay.io/cellgeni/starsolo:v4.1` |
| `local/reprocess10x/collate` | REPROCESS10X_COLLATE | `quay.io/biocontainers/python:3.12.12` |
| `local/reprocess10x/runreport` | REPROCESS10X_RUNREPORT | `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2` |

`cellgeni/starsoloqc` and `local/reprocess10x/{sra2fastq,starsolo,loaddata}` are not included
anywhere. In particular the wired-up SRA dump is `cellgeni/sra2fastq`, not the local module of the
same name; check the `include` lines before editing one:

```bash
grep -rn "include {" subworkflows/ workflow/ main.nf
```

Module scripts in `resources/usr/bin/` are on the task's `PATH` (`nextflow.enable.moduleBinaries`).
To publish a new version of a registry module, bump `version` in its `meta.yml`, then:

```bash
nextflow module validate modules/cellgeni/wget10x
nextflow module publish modules/cellgeni/wget10x -dry-run
nextflow module publish modules/cellgeni/wget10x
```

{: .important }
Editing a module changes its task hash, so the next `-resume` re-runs every task of that process
and everything downstream. Adding a key to `meta` does the same to whichever process consumes it.
A `label` or a resource directive does not.

## Process configuration

Resources follow nf-core ([Resources and queues](guides/resources.md)):

- Give a new process a `label` (`process_single`, `process_high`, `process_long`; add
  `process_network` if it fetches from the internet).
- Put only what the label cannot say in `configs/<process>.config`, wired in with an
  `includeConfig` line in `nextflow.config` after `base.config`. Each `withName` pattern starts
  with `.*`, so it still matches inside a subworkflow (`REPROCESS10X:DOWNLOAD10X:WGET10X`).
- **Never set a `queue`** outside `configs/sanger.config`: the site profile chooses it from the
  resolved memory and time, so every process needs a realistic `time`.
- A config closure cannot see `lib/` classes; compute in the workflow, put the result in `meta`,
  and read it in the directive — as `meta.resources` does for STARsolo.

The error strategy in `nextflow.config` retries every task once whatever its exit code. Its
`task.attempt == 1` clause is deliberate and load-bearing; see [Retries](guides/resources.md#retries).

## Parameters and validation

To add or change a parameter:

1. Add its default to `params {}` in `nextflow.config`. Use `null` for "not set".
2. Describe it in `nextflow_schema.json`: type, description, `help_text`, allowed values. The
   schema is tracked only through a `!/nextflow_schema.json` exception to the `*.json` rule in
   `.gitignore`.
3. A rule across parameters or profiles goes in `main.nf`, after `validateParameters()`, beside
   the existing checks — and in `MAIN_NF_RULES` in `docs/_scripts/schema_to_md.py`.
4. Regenerate the [parameter reference](reference/parameters.md):

   ```bash
   python3 docs/_scripts/schema_to_md.py
   ```

## Tests

The suite uses [nf-test](https://www.nf-test.com) and covers FETCH10XMETA against live GEO, SRA,
ENA and BioStudies metadata, so it needs outbound network access. Submit it to LSF:

```bash
bsub < tests/scripts/run_tests.bsub
bsub -env "all, ARGS=--tag relation-recovery" < tests/scripts/run_tests.bsub   # one tag
bsub -env "all, ARGS=--update-snapshot"       < tests/scripts/run_tests.bsub   # re-record snapshots
bsub -env "all, ARGS=--ci"                    < tests/scripts/run_tests.bsub   # fail on a snapshot mismatch
```

The job is only the driver; Nextflow submits the test processes to the `transfer` queue. Results
land in `logs/nf-test/testOutput.<jobid>.log`. Tags: `geo`, `arrayexpress`, `bioproject`, `enafq`,
`orifq`, `bam`, `sra`, `mouse`, `subset`, `relation-recovery`, `stub`. The full suite takes about
ten minutes; [the module README](https://github.com/cellgeni/nf-reprocessing-public-10x/blob/main/modules/cellgeni/fetch10xmeta/README.md#tests)
says what each test covers.

To run the tests directly:

```bash
module load cellgen/nf-test/0.9.5 cellgen/nextflow/26.04.6
cd modules
nf-test test cellgeni/fetch10xmeta/tests/main.nf.test --config ../nf-test.config --profile local
```

{: .important }
`cd modules` is not optional. nf-test walks the entire launch directory to build a dependency
graph before running anything, regardless of `testsDir`; from the repository root that walk never
returns, because `nf-work` holds hundreds of thousands of files.

`--profile local` runs the test processes on the current host; `--profile lsf` submits them to the
`transfer` queue.

## This site

The site is built by GitHub Pages from `main:/docs` with the
[Just the Docs](https://just-the-docs.com/) theme, pinned in `docs/_config.yml`; a page merged to
`dev` appears with the next release. Pages are plain Markdown with front matter for the sidebar
(`title`, `parent`, `nav_order`). Link between pages as `page.md`; link to anything outside
`docs/` with a full GitHub URL, since a relative link out of `docs/` is a 404 on the site.

{: .warning }
Do not move or rename the knowledge-base pages (`archive-pathologies.md`, `failure-modes.md` and
the rest at the root of `docs/`), and do not switch on `permalink: pretty`. Their `.html` URLs are
linked from upstream reports and post-mortems.

To preview locally with Ruby and Bundler:

```bash
cd docs
bundle install
bundle exec jekyll serve
```
