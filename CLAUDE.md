# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow DSL2 pipeline that loads and reprocesses public 10x datasets from GEO, SRA, ENA
or ArrayExpress: fetch metadata → download (FASTQ / BAM / SRA) → normalise to Cell Ranger
FASTQ naming → align with STARsolo or Cell Ranger. It runs on the Sanger farm via LSF and
Singularity; [nextflow.config](nextflow.config) hard-codes farm paths for whitelists,
references and the Singularity cache.

[README.md](README.md) has the user-facing parameter table, input format and output layout.
Most of the accumulated operational knowledge lives in two places, both worth reading before
any non-trivial change:

- [.claude/skills/reprocess-debug/](.claude/skills/reprocess-debug/) — the `reprocess-debug`
  skill: collecting a finished run's logs (`bin/collect_run_logs.sh`), classifying and
  attributing failures (`bin/triage.py`), exit-code meanings per process, how to verify claims
  against the `nf-work` dirs where they survive (batches 1-5 were deleted 2026-09-14),
  known-good regression accessions, and open issues found
  but not fixed. Also `references/run-index.md` — every run's name, session, archive path and
  failure counts, so none of that has to be hunted for. Supersedes the old `docs/agent_debug.md`,
  which stays untracked and local.
- [docs/](docs/) — the user-facing knowledge base, tracked and browsable on GitHub:
  [`archive-pathologies.md`](docs/archive-pathologies.md) (per-accession defects found in public
  submissions, and the separate list of correctly-labelled data we rejected late),
  [`reporting-upstream.md`](docs/reporting-upstream.md) (how to report one to GEO/SRA/ENA),
  [`failure-modes.md`](docs/failure-modes.md) (exit codes, baseline distribution, open problems),
  [`post-mortems.md`](docs/post-mortems.md) (every published run write-up and its link — the
  reports go out as Claude Artifacts, browser-readable, and this file is where the URLs are
  tracked; read a link from here rather than listing artifacts, and add one here when you
  publish).
- [docs/10x_chemistry_reference.md](docs/10x_chemistry_reference.md) — the chemistry/geometry
  table the inference scripts implement: which whitelist means which chemistry, what is
  `layout_only` vs unique, and which layouts (ATAC, feature-barcode, V(D)J, Flex) must be
  rejected from a GEX path rather than guessed at.

## Commands

Nextflow version matters: the manifest requires `>=26.04.1` and the farm default on PATH is
25.04.4, so any run launched from the repo root without a module load is refused.

```bash
module load cellgen/nextflow/26.04.6

# a run (the usual route is bsub — see scripts/run_reprocess.bsub for the live invocation)
nextflow run main.nf --datasets data/tables/batches/batch6.csv \
  --outdir results/batch6 --default_specie human --starsolo -ansi-log false

bsub < scripts/run_reprocess.bsub     # edit the nextflow line in the script first
```

Tests (nf-test, live network — needs a head node or the `transfer` queue):

```bash
bsub < tests/scripts/run_tests.bsub
bsub -env "all, ARGS=--tag relation-recovery" < tests/scripts/run_tests.bsub
bsub -env "all, ARGS=--update-snapshot"       < tests/scripts/run_tests.bsub

# directly, outside LSF — note the cd, it is required
module load cellgen/nf-test/0.9.5 cellgen/nextflow/26.04.6
cd modules
nf-test test cellgeni/fetch10xmeta/tests/main.nf.test --config ../nf-test.config --profile local
```

**`cd modules` is not optional.** nf-test walks the entire launch directory to build a
dependency graph before running anything, regardless of `testsDir`; from the repo root that
walk never returns because `nf-work` holds hundreds of thousands of files on Lustre.

Tags available: `geo`, `arrayexpress`, `bioproject`, `enafq`, `orifq`, `bam`, `sra`, `mouse`,
`subset`, `relation-recovery`, `stub`. The suite currently covers `FETCH10XMETA` only.

Iterating on chemistry inference is much faster by calling the Python directly against a
persisted work dir than by running the pipeline — see §9 of `docs/agent_debug.md`. The work
dirs only survive from batch 6 onward.

## Architecture

```
main.nf                     param validation, publish:/output: blocks, versions + QC collection
└── workflow/main.nf        REPROCESS10X — the real orchestrator
    ├── modules/cellgeni/fetch10xmeta          → links.tsv, one row per run
    ├── subworkflows/local/download10x         → WGET10X → {BAM2FASTQ,SRA2FASTQ} → RENAME10XRUN → RENAME10XSAMPLE
    ├── subworkflows/local/starsolo10x         → instantiated twice, once per species
    ├── modules/cellgeni/cellranger/count      → likewise, human + mouse
    └── modules/local/reprocess10x/mappingqc   → one QC job per dataset, across all its species
```

`modules/cellgeni/*` are registry modules (imported bare, e.g. `from 'cellgeni/fetch10xmeta'`,
and carrying a `.module-info` checksum); `modules/local/*` are repo-local and imported by
relative path. Per-process resources live one-file-per-process in [configs/](configs/), each
wired in by an `includeConfig` line at the bottom of `nextflow.config`.

### The load-bearing ideas

**`links.tsv` is read whole, not row by row.** In
[subworkflows/local/download10x/main.nf](subworkflows/local/download10x/main.nf) the whole
file is `flatMap`ped so that species, per-sample run counts and the dataset's alignable-sample
count are all settled before any of them becomes a grouping key. None of the three can be
decided from a single run's row — a sample's runs routinely carry missing or conflicting
species annotations.

**`meta.sample_id` and `meta.dataset_id` are `GroupKey`s, not strings.** They carry the group
size needed by the `groupTuple` downstream, so `.getGroupTarget()` is required to read the
accession back. `main.nf` has `unwrapGroupKeys()` for exactly this, applied on every channel
before publishing. `dataset_id`'s group size is the dataset's *alignable* sample count, which
is what lets `REPROCESS10X_MAPPINGQC` fire once per dataset rather than once per species.

**Chemistry is resolved once, per sample, and pushed into `meta`.** `RENAME10XRUN` writes a
per-run `chemistry.json`; `resolve_run_chemistry` collapses a sample's reports to one id (or
`null` if the runs disagree), and two small mappers translate it into `meta.chemistry` for
Cell Ranger's `--chemistry` and `meta.wl` for `starsolo 10x --wl`. Both keys are set *only*
when inference actually decided something, which is why `alignerIndexRow()` in `main.nf`
exists — the CSV index writes its header from the first record, so a sample missing a key
would otherwise emit a short row. Read the comments in `download10x/main.nf` before changing
which chemistry ids are forwarded: the `layout_only` ids are safe for STARsolo's `--wl` but
deliberately not passed to Cell Ranger.

**Runs of mixed origin are merged before `RENAME10XSAMPLE`.** BAM-derived and FASTQ/SRA-derived
runs of one sample join in a single channel, so each sample produces exactly one renamed set
and one aligner job. Species is deliberately *not* part of that grouping key — keying on it
would split a sample into two STARsolo dirs that then collide in MAPPINGQC.

**Species branching.** Samples resolve to `Homo sapiens` / `Mus musculus` / `UNKNOWN`; only
the first two reach an aligner, the rest are logged and dropped. `--default_specie` supplies
the fallback for missing or contradictory metadata; `--no_infer_specie` forces it everywhere.

## Gotchas

**`.gitignore` swallows most of what you write.** `data`, `scripts`, `examples`, `*.csv`,
`*tsv`, `*list`, `*.json`, `*.log`, `results*`, `reports`, `*bsub` are all ignored. `git status`
staying clean after you add a file is expected, not a bug — and note the extension rules have no
leading slash, so they match at **any** depth. A `.tsv` or `.json` anywhere, including under
`.claude/`, is ignored. New tracked files have to be `.md` or `.sh`.

`docs/` used to be ignored outright, which is why `docs/10x_chemistry_reference.md` is tracked
(it predates the rule) while anything added beside it vanished. The rule is now
`docs/agent_debug*.md`, so the knowledge base can live there. `.claude/` was never ignored —
the skill was simply untracked until it was committed; only `.claude/settings.local.json` is
excluded.

**Every task retries once regardless of exit code, and that is deliberate.** The
`errorStrategy` in `nextflow.config` reads `task.attempt == 1 || (exitStatus in [...] && attempt
< 3)`; the first clause is load-bearing and must not be "simplified" away. Transient failures
are not identifiable from an exit code — in batch 6 a STARsolo task died at exit 1 on a nonsense
read length off a partly-staged FASTQ and succeeded on its second attempt, and exit 1 is not in
the retriable set. Losing a good sample costs more than rerunning a doomed one. The retry also
means a task is only ignored after failing *twice*, so the first attempt's work dir and stderr
survive for triage.

The consequence for debugging: **`attempt=2` is not diagnostic of anything.** It is a property
of the config, not a signal about the failure — 1795 of 1814 August 2026 failures sat there.
What matters is whether any attempt later succeeded.

**`errorStrategy 'ignore'` means failures never reach the pipeline exit status**, and there is
no failure manifest, so post-mortems have to be mined out of the logs.

**Two `sra2fastq` modules exist and only one is wired up.** `download10x` imports
`modules/cellgeni/sra2fastq`; `modules/local/reprocess10x/sra2fastq/` — the more careful
implementation — is dead code. Check the `include {` lines before editing any module:

```bash
grep -rn "include {" subworkflows/ workflow/ main.nf
```

**Input tables must be tab-separated whatever the extension.** `data/tables/batches/*.csv` are
TSV; `workflow/main.nf` parses with `splitCsv(sep: '\t')` and `sample_id` itself holds
comma-separated accessions.

**Editing a process invalidates its task hash**, so touching an aligner module re-runs every
aligner task on the next `-resume`. Adding a key to `meta` does the same to whichever process
consumes it.

**`REPROCESS10X_BAM2FASTQ`'s stub is broken** — it touches `<id>.fastq.gz` but declares
`path("fastqs/*")`, which blocks stub-testing any BAM-origin path.

## Farm paths

Whitelists `/nfs/cellgeni/STAR/whitelists`, STAR references
`/nfs/cellgeni/STAR/{human,mouse}/2020A/index`, Cell Ranger references
`/software/cellgen/cellgeni/refdata_10x/`, Singularity cache
`/nfs/cellgeni/singularity/images/`. Use these directly — `find /` will exceed the tool
timeout.

Finished runs are archived to `/nfs/cellgeni/reprocessing-runs/batch<N>/<run_name>/` by
`bin/archive_run.sh` — tool stderr, manifests, trace, Nextflow reports, LSF driver logs, each
with a `manifest.txt`. Scratch is not backed up; this is the only durable copy.
`/nfs/cellgeni/projects/` is **not** writable by ordinary group members (owned by `cellgeni-su`,
mode 755), which is why the archive sits at the top level instead.
