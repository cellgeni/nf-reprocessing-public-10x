# Where the failure data lives

Distilled from `docs/agent_debug.md` §2 and §11, re-verified against disk in September 2026.
Paths are relative to the repo root; every command below assumes you are there.

## §history — resolving a run

`.nextflow/history` is the authoritative map from batch to run. Tab-separated, **no header**,
7 columns:

```
timestamp(yyyy-MM-dd HH:mm:ss) | duration | run_name | status(OK/ERR) | script_id | session_uuid | full command line
```

Column 7 holds the `--datasets data/tables/batches/batchN.csv` path, which is the only link
between a batch number and a run name — `scripts/run_reprocess.bsub` does **not** pass `-name`,
so every run gets a generated `adjective_surname`.

```bash
awk -F'\t' 'index($7,"batch6.csv"){printf "%s  %-22s %-4s %s\n",$1,$3,$4,$2}' .nextflow/history
```

**A `-resume` keeps the session uuid and takes a new run name.** Batch 5 is
`awesome_neumann` (ERR, 2026-09-07) then `tender_brattain` (OK, 2026-09-08) on the same
session `f97ddf5e…`. Triaging the wrong one of a pair shows failures that were later fixed by
the resume, or hides them. `bin/collect_run_logs.sh` refuses to guess and makes you pass
`--run` or `--latest`.

The LSF driver job id is a third identity: `logs/lsf/reprocessOutput<JOBID>.log`. Match it to a
run by the `Launching \`main.nf\` [<run_name>]` line on line 2 of that file.

## §files — per-run artifacts

| Path | What it is |
|---|---|
| `data/tables/runlogs<N>.tsv` | `nextflow log -f …` export, 36 columns, header included |
| `data/tables/failedjobs<N>.tsv` | one row per failed task, last attempt: `attempt exit name tag hash workdir` |
| `data/tables/failed<N>.log` | per-task stderr in `### <tag> (<workdir>) ###` blocks |
| `data/tables/failures<N>.tsv` | manifest: `process tag exit attempts_failed recovered workdir first_error` |
| `data/tables/batches/batch<N>.csv` | the input table for that batch (TSV despite the extension) |
| `data/tables/upload_batch<N>.csv` | delivery manifest, genuinely comma-separated: `id,dataset_id,path` |
| `results/batch<N>/` | outdir |
| `logs/lsf/reprocess{Output,Error}<JOBID>.log` | driver stdout/stderr from the bsub |
| `reports/execution_trace_<ts>.txt` | live trace, **16 columns** — a different, smaller schema than `runlogs<N>.tsv` |
| `reports/failures_<ts>.tsv` | one hand-made manifest from batch 3; same schema as `failures<N>.tsv` |
| `nf-work/<hh>/<hash…>/` | task work dirs — **these persist**, see `verify.md` |

`<N>` is the batch number, no zero-padding, appended straight to the stem. **Batch 1 is the
exception**: its trio is unsuffixed and lives in `data/`, not `data/tables/`
(`data/failedjobs.tsv`, `data/failed.log`, `data/runlogs.tsv`, `data/searchlist.tsv`).

What actually exists is patchy — as of September 2026, `{failedjobs,failed,runlogs}{2,3,5}` plus
the batch-1 trio. Batches 4 and 6 finished with no triage artifacts at all. Check before
assuming; `docs/agent_debug.md` listed a set that was already stale.

## §reference — cross-run tables

| Path | What it is |
|---|---|
| `data/tables/allhumandatasets.tsv` | `dataset_id \t sample_id(,-sep)` — 4602 datasets, 43344 unique samples |
| `data/tables/datasets.tsv` | same schema, only 266 datasets — a subset, never use it alone |
| `data/tables/searchlist.tsv` | **no header**, 5 cols: `run_id, species, urls(;-sep), type, sample_id` |
| `data/tables/All_10x.sample_table.tsv` | 10 MB master sample table |
| `data/tables/datasets.csv` | **not** a sample list — an iRODS archive manifest (`type,path,size,checksum`) |

`searchlist.tsv` is per-batch and gets overwritten, so an old one silently mis-attributes a new
run's runs. `bin/triage.py` prints which tables it used and how many tasks stayed unattributed;
a large unattributed count usually means the searchlist belongs to a different batch.

Processed samples in the iRODS manifest are the depth-6 collection paths:

```bash
awk -F',' '$1=="collection"{print $2}' data/tables/datasets.csv \
  | awk -F'/' 'NF==6{print $6}' | sort -u
```

## §sizes — what will kill your session

* `data/tables/failed2.log` is **793 MB**. `data/failed.log` is 12 MB / 227k lines for 1814
  jobs, almost all of it embedded `.command.run` boilerplate. Never `cat`, never grep unfiltered.
* `.nextflow.log` is 35 MB, with nine rotated siblings.
* `nf-work/` holds 259 hash-prefix dirs and hundreds of thousands of files on Lustre; `.lineage/`
  is comparable. `find` across either does not return inside the tool timeout, and neither does
  `find /`.
* The Bash tool timeout is 120 s. Anything scanning every work dir, or decompressing a multi-GB
  FASTQ, must be launched with `run_in_background: true` (or `nohup … &`) and polled. Do **not**
  wrap it in `timeout 900 …` in the foreground — the tool kills the call first.

## §schemas

`runlogs<N>.tsv`, 36 columns, alphabetical:

```
accelerator accelerator_type attempt complete container cpu_model cpus disk duration exit hash
hostname memory module name native_id pcpu peak_rss peak_vmem pmem process queue read_bytes
realtime rss scratch start status submit tag task_id time vmem wchar workdir write_bytes
```

`reports/execution_trace_<ts>.txt`, 16 columns, order fixed by `nextflow.config`:

```
task_id hash process tag name status exit attempt submit duration realtime queue cpus memory
peak_rss workdir
```

`status` in either is one of `COMPLETED`, `CACHED`, `FAILED`, `ABORTED`. A resumed run's export
is mostly `CACHED`; batch 5's `tender_brattain` was 7883 CACHED, 199 FAILED, 25 COMPLETED across
8107 rows.

`data/failed.log` blocks, as produced before `bin/collect_run_logs.sh` existed:

```
### <TAG> (<workdir>) ###
<the task's own stderr — the part you want>

------------------------------------------------------------
Sender: LSF System <lsfadmin@node-…>
… entire .command.run, ~200 lines of nxf_ bash boilerplate …
Resource usage summary:
    Max Memory : …
```

Everything from `Sender: LSF System` on is noise except the `TERM_*` kill reason.
`bin/collect_run_logs.sh` strips it and keeps `TERM_*` as a single `[lsf] TERM_MEMLIMIT` line,
so the new `failed<N>.log` files are a fraction of the size and greppable.

## §batches — the input tables

`scripts/split_batches.py` splits `dataset_id/sample_id` tables into batch files:

* Datasets are never split across batches; a batch closes when the next dataset would exceed
  the cap.
* The cap counts **sample entries, not unique samples** — some GSMs appear under several GSEs
  (48823 entries vs 43344 unique) and each dataset/sample pair is a separate unit of work.
* Input order (sorted by `dataset_id`) is preserved.

```bash
python3 scripts/split_batches.py --input data/tables/allhumandatasets.tsv \
  --exclude ../nf-reprocessing-public-10x/data/tables/batch1.tsv \
  --cap 1000 --outdir data/tables/batches --suffix .csv --dry-run
```

**Content must be tab-separated whatever the extension.** `sample_id` holds comma-separated
accessions and `workflow/main.nf` parses with `splitCsv(sep: '\t')`. Files named `.csv` are
fine; comma-*delimited* files are not. CRLF is harmless — `splitCsv` strips `\r`.

Current state: `batch1.csv`…`batch49.csv`, 4341 datasets, 46791 sample entries, 41471 unique
samples, with batch1.tsv excluded.
