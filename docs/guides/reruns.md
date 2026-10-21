---
title: Reruns and cleanup
parent: Guides
nav_order: 5
---

# Reruns and cleanup
{: .no_toc }

1. TOC
{:toc}

## Failures do not stop the run

A failing task is retried (see [Retries](resources.md#retries)) and then **ignored**: the rest of
the run carries on, and the run ends with exit status 0 however many samples it lost. Two places
say what was lost:

- the run report's **Samples without a mapped output** table, with the stage, process, exit
  class and work directory of each sample's last failure ([Run report](../reference/run-report.md));
- `<outdir>/cleanup/failed_tasks.tsv` and `rerun.tsv`, with `-profile lineage_keep` (below).

[Failure modes](../failure-modes.md) says what the common exit codes mean and which failures are
the data's rather than the pipeline's.

## Resuming

`-resume` with the same work directory and launch directory brings back every task whose inputs
and code are unchanged from the cache:

```bash
nextflow run main.nf -profile sanger --datasets datasets.tsv --starsolo -resume
```

Pass the profile again — profiles are not remembered between runs, and without one the pipeline
refuses to start. To resume a particular earlier run rather than the last one, give its session
id (column 6 of `.nextflow/history`): `-resume <session-id>`.

{: .important }
A task re-runs if anything that goes into its hash changed: its script, its inputs, its
container, or a key in the `meta` map it receives. Editing a module re-runs every task of that
process on the next `-resume`, and everything downstream of it. Changing a resource request or a
label does not.

## Cleaning work directories as you go: `lineage_keep`

Raw data dominates the disk, and on the farm full batches have filled the group's 142 TB quota
more than once. `-profile lineage_keep` loads the [nf-lineage-keep](https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/plugins/nf-lineage-keep)
plugin, which empties a task's work directory once its whole downstream lineage has succeeded,
and keeps everything upstream of a failure so that it can be debugged and resumed:

```bash
nextflow run main.nf -profile sanger,lineage_keep --datasets datasets.tsv --starsolo --no_raw
```

- A cleaned work directory keeps only `.command.*` and `.exitcode`.
- Every task upstream of a final failure, and every task of the failed sample, is kept intact.
- STAR's `_STARtmp` is removed whatever the outcome; it is never needed.
- `<outdir>/cleanup/` gets `deleted.tsv` (every directory cleaned, with the bytes freed),
  `failed_tasks.tsv`, and **`rerun.tsv`** — a `--datasets` table of just the failed samples.

Re-run the failures with `rerun.tsv` and `-resume` on the same work directory: everything
upstream of each failure was kept, so those tasks come back from the cache.

The plugin is not in the Nextflow registry. Install it once with `make install` in
`plugins/nf-lineage-keep` (on the farm, from a job on the `transfer` queue); see its README.

## `--no_raw`

Publishing hard-links outputs from the work directory into `--outdir`, so emptying a work
directory frees only what was never published. `--no_raw` publishes nothing under `raw/`
(downloads, per-run and renamed FASTQs) and none of their index CSVs. Combined with
`lineage_keep`, that is what makes the FASTQs actually go once a sample is done. Without
`lineage_keep`, the files still live in `nf-work/`, and only there.
