# nf-lineage-keep

A Nextflow plugin that cleans work directories during a run without throwing away what a
failure needs. A task's work directory is emptied once its whole downstream lineage has
succeeded. Every task upstream of a final failure is kept, so the failure can be debugged
and a `-resume` of the failed items finds them all in the cache. At the end it writes a
table of the failed items in the pipeline's own input format.

It exists because none of the alternatives keep inputs of ignored failures:

| | when | failed task's inputs |
|---|---|---|
| `cleanup = true` | end of a successful run | deleted, with the cache |
| nf-boost `boost.cleanup` | during the run | deleted once the failure is ignored |
| nf-gc | during the run | kept until the consumer process closes |
| **nf-lineage-keep** | during the run | kept, with every ancestor |

Lineage (`lid://`) records were considered and rejected as the engine: they point upstream
only, appear after a task finishes, and `lineage find` walks the whole store per query.

## How it decides

The plugin is a `TraceObserverV2`. It builds the task graph from files: a task whose input
lives in another task's work directory consumes that task. From the DAG at flow start it
knows, for every process output, which processes read it and whether it can be published.

A task is **OK** when
1. it succeeded or was cached;
2. every output that may be published has been;
3. no reader of its outputs is still to come: each reading process has terminated, or it is
   listed in `singleUseProcesses` and one of its tasks already consumed this task;
4. every consumer task is itself OK — recursively, but not through `barrierProcesses`.

An OK task's work directories (all attempts) are emptied except `.command.*`/`.exitcode` and
the outputs a reader may still need — in practice, what the end-of-run barriers collect.
Those are **held** and removed one by one as their readers terminate.

A task that fails for good (its error strategy ended in `ignore`, or the run aborted) is
**pinned**, with every ancestor and every task sharing its `(dataset, sample)` key. A
sample that went ahead on partial data after one of its runs failed is therefore kept
whole. Pinned tasks are never touched. Barrier failures pin nothing.

`scratch` globs are removed as soon as any attempt ends, whatever the outcome.

Failed attempts hold nothing a reader or a `-resume` uses: a failed attempt emits no
outputs, and a `-resume` re-runs a failed task in a new directory. Two switches empty them
as soon as they end, without waiting for the lineage:

- `cleanRetriedAttempts` — an attempt that failed and will be retried;
- `cleanFailedTasks` — every attempt of a task that failed for good. Pinning is unchanged:
  the task's ancestors and its sample stay whole, so the re-run still finds them cached.

Both keep the logs, symlinks and regular files under `failedKeepBelow` (default 1 MB), which
is what a failure is debugged from. Both keep the **directory itself**. On `-resume`,
Nextflow walks a task's attempt hashes in order and only moves on to attempt n+1 while
attempt n's directory exists. If a failed attempt's directory were removed, a task whose
retry succeeded would run again, and so would everything downstream of it.

## Outputs (`<outputDir>/cleanup/` by default)

| file | content |
|---|---|
| `deleted.tsv` | one row per clean: action, process, task, hash, dir, files, `bytes` (apparent), `freed_bytes` (files with a single link). Actions: `delete` (an OK task), `delete-held` (its held outputs), `scratch`, `retried` and `failed` (failed attempts), each `dry-run-`/`dry-run` in a dry run |
| `processes.tsv` | per process output: reading processes, publishable — check this first on a new pipeline |
| `failed_tasks.tsv` | final failures with their re-run key and the work dirs kept upstream (with `cleanFailedTasks`, the failure's own dir holds only logs and small files) |
| `rerun.tsv` | `dataset_id<TAB>sample_id`, samples comma-joined, empty = whole dataset. Only written when the run completed |

With `publish mode 'link'`, deleting a work file whose published twin exists frees nothing:
`freed_bytes` says what was actually released.

## Config

```groovy
plugins { id 'nf-lineage-keep@0.2.0' }

lineageKeep {
    enabled            = true
    dryRun             = true      // default: log only. Leaving it out does NOT turn it off
    keepLogs           = true      // default
    keepFailedAttempts = false     // default: superseded failed attempts go with their task
    cleanRetriedAttempts = false   // default: true empties a retried attempt as it ends
    cleanFailedTasks   = false     // default: true empties a final failure's own attempts
    failedKeepBelow    = '1 MB'    // default: what those two leave, besides logs and symlinks
    barrierProcesses   = ['REPORT', 'COLLATE']
    singleUseProcesses = ['ALIGN', 'RENAME*']
    scratch            = ['ALIGN:{_STARtmp,**/_STARtmp}']
    metaVariable       = 'meta'     // re-run key: meta[datasetKey], meta[sampleKey], meta[idKey]
    datasetKey         = 'dataset_id'
    sampleKey          = 'sample_id'
    idKey              = 'id'
    sampleListVariable = 'sample_ids'
    outputDir          = null       // default <outputDir>/cleanup
    interval           = '30s'
}
```

Process patterns match the simple name (after the last `:`); `*` is a wildcard. Processes
matching one pattern form one group, so `STARSOLO*` treats both species aliases as one
single-use reader.

Re-run key: `dataset = meta.dataset_id ?: meta.id`; `sample = meta.sample_id`, or `meta.id`
when `dataset_id` is present without `sample_id`. A failure with no sample re-runs its
`sample_ids` subset, else the samples upstream of it, else the whole dataset.

Getting `singleUseProcesses` wrong is the one way to delete too early: a second task of the
process would find its input gone. The plugin warns when a pin reaches a task it already
cleaned.

## Build and install

JDK 21 (`module load HGI/common/jdk/21.0.12.1` on the farm; the Gradle plugin asks for a 21 toolchain). The build downloads Gradle and Nextflow, so on the farm run it on the `transfer`
queue:

```bash
make test      # Spock tests
make install   # unpacks into ${NXF_PLUGINS_DIR:-~/.nextflow/plugins}/nf-lineage-keep-0.2.0
```

## Limits

- Only the samples a run touched are in `rerun.tsv`; a sample the metadata step silently
  dropped never reached a task. Cross-check against the input table.
- A `-resume` with the *original* table re-runs every cleaned task.
- Files made by operators rather than tasks (`collectFile`) are invisible to the graph.
- Publishing is tracked for the workflow `output {}` block only. `publishDir` copies are not
  waited for, so a pipeline that publishes with `publishDir` must not use this plugin.

## Validation

`validation/` is the pipeline's shape in miniature, with no network: per-dataset metadata,
per-run download and rename, per-sample merge and align, per-dataset QC, two barriers, every
stage published by hard link, and one sample whose align fails for good. `validation/run.sh`
runs it, checks what was kept and what was cleaned, then re-runs `rerun.tsv` with `-resume`
and checks that the failed sample's upstream tasks come back `CACHED`.
