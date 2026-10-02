---
name: reprocess-debug
description: >
  Use when debugging a run of the nf-reprocessing-public-10x pipeline: collecting
  a finished run's logs, finding out what failed and why, classifying and
  attributing failures to datasets, screening a finished batch's outputs for
  samples that aligned to nothing although no task failed, checking a
  chemistry-inference or SRA2FASTQ claim against the persisted work dirs,
  writing the post-mortem report, drafting a dated issue and suggested diff in
  ./issues/ for a confirmed defect, or notifying the team by email/Slack with the
  batch summary and project progress.
  Triggers on "triage/debug batch N", "what failed in the last run", a pasted
  Nextflow run name, an nf-work path, a "terminated with an error exit status"
  line, "did anything come out wrong", a sample with cells but almost no
  features, a request for a failure post-mortem, "draft an issue"/"file this",
  or "send an update"/"email me"/"notify".
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

Five rules override everything below.

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
   work dirs for batches 1-5 were deleted (confirmed 2026-09-14), and those for batches 6-19 by
   2026-10-02 (0 of 818 failed-task dirs left)**, so walk-back works only on runs still on
   scratch; for the rest the archive under `/nfs/cellgeni/reprocessing-runs/` is the whole of the
   surviving evidence. Check the dir exists before promising to verify anything.
4. **Mind the sizes.** The Bash tool timeout is 120 s; `data/tables/failed2.log` is 793 MB;
   `.nextflow.log` is 35 MB; `find` across `nf-work/` or `/` does not return. Long collections
   go in the background and get polled. Details in `data-map.md §sizes`.
5. **A task that completed is not an output that is good, and a sample that never ran leaves no
   trace at all.** Batch 9 published 38 near-empty matrices, batch 10 published 75 and batch 11
   published 144 — 4.8%, 10.0% and **17.1%** of everything they aligned — non-GEX libraries
   accepted as GEX, with not one task failing. Batch 11 had the cleanest failure list of any
   batch (8 permanent, no dataset lost) and the worst outputs, so the two are not correlated. Batch 10 also lost 6 samples inside `FETCH10XMETA`, which exited 0, and pooled
   four experimental conditions into one published matrix that looks entirely healthy. The
   failure list shows none of this and `triage.py` never sees it. Screen the successes (§4)
   before concluding anything about a batch, and do it **while `results/<outdir>/` still
   exists** — it is deleted after upload, and with it the only record of whether the outputs
   were any good.

## 1. Route

| Stage | Read | When |
|---|---|---|
| **collect** | `data-map.md` | a run finished and has no `failures<N>.tsv` yet |
| **classify** | `classify.md` | you have the artifacts and want the shape of the failure set |
| **screen** | `verify.md §qc` | you have the failure picture and now need to know whether what *succeeded* is any good |
| **verify** | `verify.md` | before believing any claim about read structure, chemistry or a guard |
| **report** | `report.md` | the answer is going to a person, not just the terminal |
| **file** | `issues.md` | a finding is a confirmed, actionable defect — draft it for the tracker |
| **notify** | `notify.md` | a debug session is finished — save the LSF log, refresh the progress counter, email/Slack the team (`--team`) |
| **compare** | `known-issues.md` | is this run's profile normal, and is this bug already known |
| **look up** | `run-index.md` | which run was batch N, where is its evidence archived, is there a post-mortem — and `docs/post-mortems.md` for its link |

| Signal | Go to |
|---|---|
| "where is the batch-N post-mortem", "send me the link" | `docs/post-mortems.md` — do not go listing artifacts |
| "which run was batch N", "where is its evidence" | `run-index.md` |
| "triage batch N", "what failed in the last run" | §2, §3, then §4 |
| "did anything come out wrong", a batch with few or no failures | `verify.md §qc` |
| a sample with a plausible cell count but almost no features | `verify.md §qc` |
| a run name (`spontaneous_ampere`) or an LSF job id | `data-map.md §history` |
| an `nf-work/<hh>/<hash>` path | `verify.md §walkback` |
| `terminated with an error exit status (N)` | `classify.md §exits` |
| a run still in flight | `classify.md §console` |
| "why was this sample rejected", a chemistry argument | `verify.md §whitelist`, then `docs/10x_chemistry_reference.md` |
| "write it up", "post-mortem", "share the findings" | `report.md` |
| "file this", "draft an issue", "open a ticket", a confirmed defect worth fixing | `issues.md` |
| "email me", "notify ab76", "send an update", "post to slack" | `notify.md` |
| "who gets the emails", "X isn't getting the reports" | `notify.md` — recipients are `notify.py`'s `TEAM` |

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

## 4. Screen what succeeded

The failure list says nothing about the samples that aligned, nor about the ones that never
reached a task. Each of the last three batches had its largest finding here: 38 samples in
batch 9, **75 in batch 10** and **144 in batch 11** — hashing, antibody, CMO, guide-RNA and
V(D)J libraries accepted as GEX, plus an entire 97-sample bulk amplicon dataset that is not 10x
at all — plus, in batch 10, 6 samples dropped by the metadata stage without a word. None of it
appears in any failure list. Batch 11 is the case to remember: 8 permanent failures, no dataset
lost, and 17.1% of everything it published was not gene expression.

```bash
# what aligned to nothing. De-duplicate on Sample first — duplicate dataset rows inflate this.
awk -F'\t' 'NR>1 && !seen[$2]++ && $17+0 < 0.05 && $8+0 < 100 {print $1, $2, $14, $17, $8}' \
    results/batch10/mapping_qc_stats.tsv

# what never got as far as a task: requested samples missing from the run's own links.tsv
python3 - <<'EOF'
import csv, os
B = "10"
# batch22+ live in batches_deduplicated/; the never-run batches/batch22-49 share the names
bt = next(p for p in (f"data/tables/batches_deduplicated/batch{B}.csv",
                      f"data/tables/batches/batch{B}.csv") if os.path.exists(p))
for r in csv.DictReader(open(bt), delimiter="\t"):
    ds  = r["dataset_id"]
    req = {s.strip() for s in r["sample_id"].split(",") if s.strip()}
    p   = f"results/batch{B}/metadata/{ds}/links.tsv"
    if not os.path.exists(p):
        print(f"{ds}: no links.tsv at all ({len(req)} samples)"); continue
    seen = {l.split("\t")[4].strip() for l in open(p) if len(l.split("\t")) >= 5}
    if req - seen:
        print(f"{ds}: {len(req & seen)}/{len(req)} emitted, missing {sorted(req - seen)}")
EOF
```

**Both terms of the first screen are load-bearing.** `exon_u < 0.05` alone called a
single-nucleus sample junk in batch 10 — `GSM6729514`, `exon_u` 0.043 but 583 median features
and `full_u` 0.47. `Med_nFeature < 100` removes exactly that and nothing else. Do **not** reach
for `full_u` as the second term: it would clear 7 of batch 9's confirmed hashing libraries, which
reach `full_u` 0.38 on 1-5 features. Both numbers are screening values from three batches, not
verdicts — read `verify.md §qc` before acting on a hit.

**The first screen is precise and not sensitive — run the title scan too.** It has never
produced a false positive, but it only finds libraries that map to *nothing*, and that is not
the whole population. In batch 11 it flagged 42 of a measured 144 non-GEX samples. It missed
23 V(D)J (TCR) libraries that map at 73-90% with up to 246 median features, and 79 of the 97
samples of a bulk amplicon dataset that maps exonically at 0.05-0.24. Do not read a small
first-screen result as a clean batch.

```bash
# The libraries GEO names in the title. Zero false positives in batch 11 (47 of 844),
# where the first screen found 42 of 144. Token-delimited on purpose: a substring match
# also pulls in Ascites_IG_1, "scRNA CRISPRa ... S1" and "spike-specific B cells", all
# healthy GEX.
for f in results/batch11/metadata/*/*_family.soft; do
  awk -v d="$(basename "$(dirname "$f")")" \
      '/^\^SAMPLE = /{g=$3} /^!Sample_title/{sub(/^!Sample_title = /,""); print d"\t"g"\t"$0}' "$f"
done | grep -Ei '(^|[[:space:]_-])(TCR|BCR|VDJ|CITE|HTO|ADT|CSP|CMO|hashtag|hashing|gRNA|sgRNA)([[:space:]_-]|$)|	BC[0-9]+$'
```

Confirm a hit before calling it: the matrix itself is decisive. A V(D)J library's top genes are
`TRBV`/`TRAV` (or `IGHV`) segments — count UMIs per gene in
`starsolo/<GSE>/<GSM>/output/Gene/filtered` rather than arguing from the title alone.

`Med_nFeature < 100` on its own catches 131 of batch 11's 144 with one false positive, but it is
**not** a safe replacement for the conjunction: on batch 10 it adds 28 samples of which at least
27 are real shallow libraries (GSE219098 ×22 at `exon_u` 0.53). Use it as a review flag.

**The second screen catches a different failure entirely**, and it is the one with no other
symptom: batch 10's `FETCH10XMETA` emitted 2 of GSE218936's 8 requested samples, exited 0, and
aligned the survivor against all four runs of its group — publishing four experimental conditions
pooled under one sample's name, at 96% mapping and 1865 median features. Nothing about that
output looks wrong. Batch 9 is clean on this check; batches 1-8 cannot be checked.

Two operational traps:

- **De-duplicate before counting.** `mapping_qc_stats.tsv` carries one row per dataset unit, not
  per sample, so duplicate dataset accessions double-count — batch 9, 917 rows for 784 samples;
  batch 10, 973 for 748. Batch 10's 75 flagged samples first appear as 97 rows.
- **Do this before the batch is cleaned up.** `results/<outdir>/` is deleted after upload; only
  batches 9 and 10 still had one in September 2026. `archive_run.sh` (§10) preserves
  `mapping_qc_stats.tsv` and the run's `links.tsv`, but only if it runs first.

## 5. Verify before you conclude

Anything you are about to state as fact about read geometry, barcode match rate, chemistry, or
which stage really caused a failure: measure it. `verify.md` has the whitelist probe, the
staged-symlink walk-back, how to run the inference script directly against a work dir, and the
regression accessions to re-check after touching inference.

## 6. Report

If the answer is going to a person, it goes out as an Artifact, not as scrollback. `report.md`
has the procedure, the house palette and layout, the three-state theme trap, and the
reconciliation checks to run before publishing. Load the `artifact-design` skill first.

Relay the headline findings in the chat reply as well as linking the page.

## 7. Draft an issue for every confirmed defect

Not conditional on writing a report — a session with no post-mortem still drafts one of these
for any finding that is a real, actionable bug rather than a correct rejection. `issues.md` has
the trigger condition, the naming (dated, since sessions stack up between reviews), the issue
template, and the safe procedure for generating a real `git diff` without leaving the working
tree modified.

**Never run `gh issue create` or apply the suggested diff yourself.** Both files are drafts for
the user to review; say so plainly in the chat reply, with the paths.

## 8. Save the LSF log and notify

Closing step for every debug session, whether or not it produced a report or an issue draft.
`notify.md` has the three parts:

```bash
.claude/skills/reprocess-debug/bin/save_lsf_log.sh --batch 7          # onto NFS now, don't wait for §10
.claude/skills/reprocess-debug/bin/track_progress.py --batch 7 --run elegant_lamarr --json /tmp/progress.json
.claude/skills/reprocess-debug/bin/notify.py --team --subject "..." --body-file /tmp/body.txt --slack
```

Email always sends (`smtplib.SMTP("localhost")`, no credentials needed). **Use `--team`, not a
hand-typed `--to`** — it expands to `notify.py`'s `TEAM` list, the one place the recipients are
defined; a hand-typed address is how ap41 silently got nothing for batches 6-9. Slack only sends if
`SLACK_BOT_TOKEN`/`CHANNEL_ID` resolve to something real — **never go looking for that token
yourself or fabricate one**; if it's not configured, say so and send the email anyway. Its
`slack_sdk` dependency is borrowed per-call through `uvx`, because this python has user
site-packages disabled and `pip install --user` therefore installs something unimportable — see
`notify.md`. Report what `notify.py`'s stderr actually says; "posted to channel" is the only line
that means it went out.

**Read `notify.md §2` about the progress percentage before quoting it anywhere.** Since
2026-10-02 its denominator is the deduplicated to-do list, `data/tables/batches_deduplicated/`
(batch22-35, 24,864 samples), not the March 2026 project-wide target. It started at 0% and
counts how much of what was left is done. It is not comparable with rows before 2026-10-02,
which were measured against 44,494 samples.

## 9. Files

| Path | What |
|---|---|
| `bin/collect_run_logs.sh` | run/batch → `runlogs`, `failedjobs`, `failed.log`, `failures` manifest |
| `bin/triage.py` | manifest → classified, attributed, verdicts, `--json <path>` |
| `bin/archive_run.sh` | a finished run → a complete record under `/nfs/cellgeni/reprocessing-runs/` |
| `bin/resolve_run.sh` | sourced helper: batch → run name, LSF job id, report stamp |
| `bin/save_lsf_log.sh` | copies a run's LSF driver log to NFS immediately, ahead of the full archive |
| `bin/track_progress.py` | progress vs. the deduplicated to-do list (`data/tables/batches_deduplicated/`) — iRODS + local results, appends to `progress-counter.tsv` |
| `bin/notify.py` | email (`--team`, or `--to` for a one-off address) + Slack (if configured) — never fabricates a missing credential; `TEAM` is the recipient list |
| `bin/notify_run_done.sh` | the unattended "run finished, N left on Lustre" Slack ping, called from `scripts/run_reprocess.bsub` — not part of a debug session |
| `references/data-map.md` | where everything lives, schemas, naming, size traps, batching |
| `references/classify.md` | rule semantics, exit codes, attribution joins, live-run triage |
| `references/verify.md` | work-dir probes, guards and thresholds, regression accessions — and `§qc`, screening a finished batch's outputs for samples that aligned to nothing |
| `references/report.md` | the Artifact post-mortem |
| `references/issues.md` | drafting a dated issue + suggested diff in `./issues/` for a confirmed defect |
| `references/notify.md` | saving the LSF log, refreshing the progress counter, emailing/Slacking the team |
| `references/known-issues.md` | baseline distributions, per-batch history, open bugs |
| `references/run-index.md` | every run: name, session, failure counts, archive path |

Related, outside this skill: `CLAUDE.md` for the pipeline architecture, and `docs/` for the
user-facing knowledge base — `10x_chemistry_reference.md` (what each whitelist means and which
layouts must be rejected), `archive-pathologies.md`, `reporting-upstream.md`,
`failure-modes.md`, and `post-mortems.md` — the index of record for every published report's
link, the one place to read a URL from or write a new one to.

## 10. Archive the run when you are done

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
