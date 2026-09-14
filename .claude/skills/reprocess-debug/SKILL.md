---
name: reprocess-debug
description: >
  Use when debugging a run of the nf-reprocessing-public-10x pipeline: collecting
  a finished run's logs, finding out what failed and why, classifying and
  attributing failures to datasets, checking a chemistry-inference or SRA2FASTQ
  claim against the persisted work dirs, or writing the post-mortem report.
  Triggers on "triage/debug batch N", "what failed in the last run", a pasted
  Nextflow run name, an nf-work path, a "terminated with an error exit status"
  line, or a request for a failure post-mortem.
---

# Debugging a reprocessing run

Router. Work out which stage you are at, then read **only** those reference files.

The pipeline itself is described in `CLAUDE.md`; this skill is about its failures. User-facing
background lives in `docs/` — [`failure-modes.md`](../../../docs/failure-modes.md) for exit
codes and open problems, [`archive-pathologies.md`](../../../docs/archive-pathologies.md) for
per-accession defects, [`reporting-upstream.md`](../../../docs/reporting-upstream.md) for
sending one to GEO/SRA/ENA.

> **Local-only paths.** `scripts/` and `docs/agent_debug*.md` are gitignored working files that
> exist only on the machine this skill was written on. This skill replaces `agent_debug.md`; the
> "Distilled from `docs/agent_debug.md` §N" lines at the top of each reference file are
> provenance notes, not links you can follow. `scripts/triage3.py`, `scripts/split_batches.py`
> and `scripts/run_reprocess.bsub` are referenced the same way — if you do not have them, the
> surrounding instructions still stand but their commands will not run.

Four rules override everything below.

1. **`attempt` is not diagnostic, and neither is a run's `OK` status.** The `errorStrategy` in
   `nextflow.config` retries every task exactly once whatever the exit code — **deliberately**,
   so that a transient failure at an unretriable exit code is not thrown away; see
   `classify.md §exits` for the evidence, and do not report it as a precedence bug. The
   consequence for triage is that 1795 of 1814 August 2026 failures sat at `attempt=2` and it
   meant nothing. And
   `errorStrategy 'ignore'` keeps failures out of the exit status entirely: batch 6 is recorded
   `OK` in `.nextflow/history` with 7 permanent failures. Never report a run as clean because
   Nextflow said so.
2. **Aggregate by dataset before reporting anything.** Failures concentrate hard — the top 5
   datasets were 82% of August 2026's 1814, one alone was 48%. A flat accession list hides that,
   and a fix aimed at the flat list aims at the wrong thing.
3. **Measure, do not infer — where you still can.** Several confident log-based conclusions were
   wrong until checked against the FASTQs, because a failure routinely surfaces two stages
   downstream of its cause wearing an unrelated message. See `verify.md §walkback`. **But the
   work dirs for batches 1-5 were deleted (confirmed 2026-09-14)**, so walk-back works only from
   batch 6 onward; for earlier runs the archive under `/nfs/cellgeni/reprocessing-runs/` is the
   whole of the surviving evidence. Check the dir exists before promising to verify anything.
4. **Mind the sizes.** The Bash tool timeout is 120 s; `data/tables/failed2.log` is 793 MB;
   `.nextflow.log` is 35 MB; `find` across `nf-work/` or `/` does not return. Long collections
   go in the background and get polled. Details in `data-map.md §sizes`.

## 1. Route

| Stage | Read | When |
|---|---|---|
| **collect** | `data-map.md` | a run finished and has no `failures<N>.tsv` yet |
| **classify** | `classify.md` | you have the artifacts and want the shape of the failure set |
| **verify** | `verify.md` | before believing any claim about read structure, chemistry or a guard |
| **report** | `report.md` | the answer is going to a person, not just the terminal |
| **compare** | `known-issues.md` | is this run's profile normal, and is this bug already known |
| **look up** | `run-index.md` | which run was batch N, where is its evidence archived, is there a post-mortem — and `docs/post-mortems.md` for its link |

| Signal | Go to |
|---|---|
| "where is the batch-N post-mortem", "send me the link" | `docs/post-mortems.md` — do not go listing artifacts |
| "which run was batch N", "where is its evidence" | `run-index.md` |
| "triage batch N", "what failed in the last run" | §2, then §3 |
| a run name (`spontaneous_ampere`) or an LSF job id | `data-map.md §history` |
| an `nf-work/<hh>/<hash>` path | `verify.md §walkback` |
| `terminated with an error exit status (N)` | `classify.md §exits` |
| a run still in flight | `classify.md §console` |
| "why was this sample rejected", a chemistry argument | `verify.md §whitelist`, then `docs/10x_chemistry_reference.md` |
| "write it up", "post-mortem", "share the findings" | `report.md` |

Not covered: fixing the aligners, the reference/whitelist layout, batching strategy beyond
`data-map.md §batches`, and anything about delivery or iRODS. Say so rather than guessing.

## 2. Collect

```bash
.claude/skills/reprocess-debug/bin/collect_run_logs.sh --batch 6
```

Resolves batch → run name via `.nextflow/history`, exports the 36-field trace, reduces it to the
last failed attempt per task, and writes four files into `data/tables/`:
`runlogs<N>.tsv`, `failedjobs<N>.tsv`, `failed<N>.log`, `failures<N>.tsv`.

**Run it in the background** — `nohup … &` or `run_in_background: true`. It reads a work dir per
failed task, and a 1800-failure run will not finish inside the tool timeout.

Three things it does that hand-collection got wrong, and that you should not undo:

- **It loads `cellgen/nextflow/26.04.6` itself.** The farm PATH default is 25.04.4, which
  rejects the `accelerator` trace field and reads a 26.x cache as an empty run — no error, three
  rows, silently useless.
- **It prefers `.command.err` and truncates `.command.log` at the LSF report**, keeping only the
  `TERM_*` kill reason as one `[lsf]` line. Batch 6's `failed6.log` is 204 lines for 16 jobs;
  the old batch-1 file was 227k lines for 1814, nearly all `.command.run` boilerplate. Grepping
  the old format returns boilerplate.
- **It refuses to guess between a batch's runs.** A `-resume` keeps the session uuid and takes a
  new run name; batch 5 is `awesome_neumann` (ERR) then `tender_brattain` (OK). Pass `--run` or
  `--latest`.

If the run predates the script, or `nf-work` has been cleaned, you are stuck with whatever
`failed<N>.log` exists — go to §3 with `--failed-log`.

## 3. Classify

```bash
.claude/skills/reprocess-debug/bin/triage.py --manifest data/tables/failures6.tsv
# older runs, no manifest:
.claude/skills/reprocess-debug/bin/triage.py \
  --failed-log data/tables/failed5.log --failedjobs data/tables/failedjobs5.tsv
```

Prints process × exit × category, category × verdict, failures by dataset, route × process, the
distinct first errors with digits masked, and every `UNCLASSIFIED` row with its work dir.
`--json` writes the same aggregates for the report step.

**Read the `UNCLASSIFIED` block in full, every time.** In August 2026 the two most valuable
findings sat in a bucket of four — a silent single-mate SRA dump, and LSF memory kills wearing
exit 130. That bucket is the point of the tool, not its residue.

**Check `permanent` against `recovered`.** The retry-once bug means transient failures are
always in the raw list; every percentage should be over the permanent set. Batch 6: 16 failed
tasks, 9 self-healed, 7 real.

Rules live in `bin/triage.py` and nowhere else. Before adding one, read the blocks it would
cover, give it a `VERDICT` entry in the same change, and re-run the batch-3 regression in
`verify.md §regression` — the rule set is first-match-wins, so a new pattern can silently steal
rows from an old one.

## 4. Verify before you conclude

Anything you are about to state as fact about read geometry, barcode match rate, chemistry, or
which stage really caused a failure: measure it. `verify.md` has the whitelist probe, the
staged-symlink walk-back, how to run the inference script directly against a work dir, and the
regression accessions to re-check after touching inference.

## 5. Report

If the answer is going to a person, it goes out as an Artifact, not as scrollback. `report.md`
has the procedure, the house palette and layout, the three-state theme trap, and the
reconciliation checks to run before publishing. Load the `artifact-design` skill first.

Relay the headline findings in the chat reply as well as linking the page.

## 6. Files

| Path | What |
|---|---|
| `bin/collect_run_logs.sh` | run/batch → `runlogs`, `failedjobs`, `failed.log`, `failures` manifest |
| `bin/triage.py` | manifest → classified, attributed, verdicts, `--json <path>` |
| `bin/archive_run.sh` | a finished run → a complete record under `/nfs/cellgeni/reprocessing-runs/` |
| `bin/resolve_run.sh` | sourced helper: batch → run name, LSF job id, report stamp |
| `references/data-map.md` | where everything lives, schemas, naming, size traps, batching |
| `references/classify.md` | rule semantics, exit codes, attribution joins, live-run triage |
| `references/verify.md` | work-dir probes, guards and thresholds, regression accessions |
| `references/report.md` | the Artifact post-mortem |
| `references/known-issues.md` | baseline distributions, per-batch history, open bugs |
| `references/run-index.md` | every run: name, session, failure counts, archive path |

Related, outside this skill: `CLAUDE.md` for the pipeline architecture, and `docs/` for the
user-facing knowledge base — `10x_chemistry_reference.md` (what each whitelist means and which
layouts must be rejected), `archive-pathologies.md`, `reporting-upstream.md`,
`failure-modes.md`, and `post-mortems.md` — the index of record for every published report's
link, the one place to read a URL from or write a new one to.

## 7. Archive the run when you are done

```bash
.claude/skills/reprocess-debug/bin/archive_run.sh --batch 6 --dry-run   # always first
.claude/skills/reprocess-debug/bin/archive_run.sh --batch 6
```

Copies the tool stderr, manifests, trace, Nextflow reports, LSF driver logs and post-mortem
source to `/nfs/cellgeni/reprocessing-runs/batch<N>/<run>/` with a `manifest.txt` that says what
the run was and what is *not* in the archive. Scratch is not backed up and `nf-work` is the only
thing standing between you and an unverifiable claim.

Nothing is skipped silently; a missing or oversized file is recorded in `manifest.txt`. Then add
the run to `references/run-index.md`.
