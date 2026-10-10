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
| `data/tables/batches_deduplicated/batch<N>.csv` | the input table for **batch 22 onward** (TSV despite the extension); see §batches |
| `data/tables/batches/batch<N>.csv` | the input table for batches 1-21. Its `batch22`-`batch49` were **never run** and are superseded: same names, older and duplicated contents |
| `data/tables/upload_batch<N>.csv` | delivery manifest, genuinely comma-separated: `id,dataset_id,path` |
| `results/batch<N>/` | outdir |
| `logs/lsf/reprocess{Output,Error}<JOBID>.log` | driver stdout/stderr from the bsub |
| `reports/execution_trace_<ts>.txt` | live trace, **16 columns** — a different, smaller schema than `runlogs<N>.tsv` |
| `reports/failures_<ts>.tsv` | one hand-made manifest from batch 3; same schema as `failures<N>.tsv` |
| `nf-work/<hh>/<hash…>/` | task work dirs — **deleted for batches 1-5** (2026-09-14) **and 6-19** (by 2026-10-02); only the newest runs keep theirs, see `verify.md` |

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
| `data/tables/All_10x.sample_table.tsv` | 10 MB master sample table, **no header**, 7 cols: `sample_id, GSE(s) (,-sep or -), project, SRS, SRX(s), SRR(s), species`. 55,529 human + 48,913 mouse rows. 7,413 human samples list two GSEs; 49 rows name 2-4 GSMs for one SRA sample in col 1. The source of `batches_deduplicated/` |
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
* The Bash tool timeout is 120 s. Do **not** wrap a long command in `timeout 900 …` in the
  foreground — the tool kills the call first.

### §headnode — what may run where

The session runs on a farm head node (`farm22-head1/2`). They are for editing files, talking to
LSF and light reads, and **Arbiter** throttles anyone who uses them for compute. On 2026-10-07
the batch-23 debug session ran, from the head node and "in the background": full-file awk
scans of 30-110 GB dumped FASTQs (read-length runs, per-lane counts, distinct-read counts, up
to 5 in parallel), `gzip -t` on ~9 GB of FASTQs, and a python pass over 218 published matrices.
Arbiter reported ~300% CPU for `awk` alone and put ab76 in `penalty1` (CPU cut to 80% of 4
cores for 30 minutes, on both head nodes). `run_in_background: true` and `nohup … &` only free
the conversation; the work still runs on the head node.

**Heavy work goes through `bin/farm_run.sh`**, which submits to LSF with `bsub -K` and returns
the job's output and exit status. Launch the wrapper itself with `run_in_background: true` —
queueing alone can outlast the tool timeout — and you are notified when the job ends:

```bash
W=.claude/skills/reprocess-debug/bin/farm_run.sh
$W --name rle-SRR17720155 -- \
  'awk "NR%4==2{l=length(\$0); if(l!=p){print NR/4, l; p=l}}" nf-work/65/c3fc80…/SRR17720155_2.fastq'
$W --name matstats --mem 8G -- python3 logs/lsf/debug/matstats.py cands.txt logs/lsf/debug/matstats.json
$W --name collect23 --mem 8G -- .claude/skills/reprocess-debug/bin/collect_run_logs.sh --batch 23
```

Defaults: queue `normal` (12 h limit), 1 CPU, 4 GB, `-W 4:00`, group `cellgeni`; outputs land
in `logs/lsf/debug/<name>.<stamp>.{out,err,lsf}`. One argument after `--` runs as a bash snippet
(pipes, globs, loops), several run as an argv. The job inherits the session's environment and
directory, so `module` and relative paths work. **It cannot see the head node's `/tmp`**, which
is where the agent's scratchpad lives: a helper script there is missing on the execution node,
and an output written there lands in that node's `/tmp` and is lost (the first LSF run of
`collect_run_logs.sh --outdir <scratchpad>` did exactly that). Keep job inputs and outputs on
Lustre, e.g. `logs/lsf/debug/`; the wrapper refuses a `/tmp` path unless `--allow-tmp`. Several independent measurements are several
background `farm_run.sh` calls, which LSF runs in parallel on separate nodes. Use
`--queue transfer` for anything that needs outbound network, `--no-wait` to submit and follow
with `bjobs -J rdebug-<name>`.

| On the head node — fine | Through `farm_run.sh` |
|---|---|
| `triage.py`, `track_progress.py` (`iquest` is permitted), `notify.py`, `archive_run.sh`, `save_lsf_log.sh` | `collect_run_logs.sh` — `nextflow log` is a JVM over the whole trace, plus one work dir per failed task |
| awk/grep over `mapping_qc_stats.tsv`, `links.tsv`, `sra.tsv`, the `failures`/`failed` files, SOFT files | anything that reads a FASTQ past its first few hundred thousand records |
| `.command.err`/`.command.log` of a handful of work dirs | `gzip -t` / `pigz -t` / `md5sum` of a whole file |
| a FASTQ's head: `zcat f.gz \| head -400000` | `wc -l`, read-length runs, per-lane or distinct-read counts over a whole file |
| one matrix's top genes | a pass over many matrices, or many work dirs |
| `git`, `bjobs`, editing | `infer_10x_run_recommended.py` against a work dir (`--sample-window 2000000` takes ~40 s CPU per run) |
|  | the `§whitelist` probe at `start=1_000_000` or on more than one file |

Rough line: **more than about a CPU-minute, more than about 1 GB read, or more than one process
in parallel → LSF.** Memory-heavy work counts too: a hash of every distinct read of a 70 M-read
file is gigabytes, so bound it (exit early past a cap, as the GSE261353 check did) and give the
job `--mem`.

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

**From batch 22 on, the input tables are `data/tables/batches_deduplicated/batch22.csv` …
`batch35.csv`**, built on 2026-10-02 by `scripts/make_dedup_batches.py`. They are the
deduplicated to-do list: every human sample in `All_10x.sample_table.tsv` that was not yet on
iRODS, not completed in local `results/`, and not in a batch already run. 24,864 samples,
2,217 datasets.

* **One entry per sample.** A sample listed under several GSEs is kept under the first listed.
  A row naming several GSMs for one SRA run keeps the first GSM. So the duplicate-dataset
  alignments measured in batches 8-21 (9-23% of the stage) should not recur, and entries equal
  unique samples.
* Samples with neither a GSE nor a project take their dataset from `allhumandatasets.tsv`
  (242 of them).
* Datasets are never split. The cap is balanced to give `round(total / 1800)` batches, here 14
  of 1,678-1,842 samples with a last one of 1,334, rather than 15 with a 198-sample remainder.
* `manifest.txt` in the folder records the sources, the iRODS cache timestamp and every
  exclusion count. It is also the denominator `bin/track_progress.py` now counts against (see
  `notify.md §2`).

```bash
python3 scripts/make_dedup_batches.py --dry-run \
  $(for b in data/tables/batches/batch{1..11}.csv data/tables/batches/batch12-13.csv \
             data/tables/batches/batch{14..21}.csv ../nf-reprocessing-public-10x/data/tables/batch1.tsv; \
    do echo --run-batch $b; done)
```

The `--run-batch` list is what was already run when it was built. Regenerating it later, with
more batches run, renumbers everything from `--start` (default 22). Write to a new `--outdir`
rather than over a folder whose batches are in flight; the script refuses an existing one.

### Batches 1-21: `data/tables/batches/`

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

State of that folder: `batch1.csv`…`batch49.csv`, 4341 datasets, 46791 sample entries, 41471
unique samples, with batch1.tsv excluded. Only 1-21 were run; 22-49 are superseded by
`batches_deduplicated/`.
