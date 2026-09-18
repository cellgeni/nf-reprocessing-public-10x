# Saving the LSF log and notifying ab76 when a batch is debugged

The last thing a debug session does, after §5 Report / §6 Issues: get the run's LSF driver log
off scratch and onto NFS, refresh the project-wide progress counter, and tell ab76@sanger.ac.uk
what happened — by email always, by Slack when it is configured. Do this for every batch you
finish debugging, whether or not a post-mortem was published.

Inspiration and reused mechanics: `/lustre/scratch124/cellgen/cellgeni/aljes/track_reprocessing`
(the `sample-tracking` CLI's `emit_reports()` — plain `smtplib.SMTP("localhost")`, no auth) and
`/lustre/scratch124/cellgen/cellgeni/aljes/reprocessing_slack_bot` (`app.py` — `slack_sdk.WebClient`
with a bot token and a channel id). Both are separate, older projects; this skill's scripts copy
their approach rather than importing them, since neither is a dependency of this repo.

## §0 Not to be confused with the run-finished ping

`bin/notify_run_done.sh` is a **different** notification with a different trigger: it is the line
after `nextflow run` in `scripts/run_reprocess.bsub`, it fires unattended the moment a run stops,
and it carries only "the run stopped, here is how it stopped, here is what is left of the
`cellgeni` Lustre quota". No triage, no post-mortem link, no progress counter — those need a
human-or-agent debug session and belong to §1-§3 below. Nothing in this file is a prerequisite
for it, and running it does not discharge the §1-§3 closing steps.

```bash
.claude/skills/reprocess-debug/bin/notify_run_done.sh --batch batch10 --exit-code "$status"
.claude/skills/reprocess-debug/bin/notify_run_done.sh --batch 10 --exit-code 0 --dry-run   # see it first
```

Slack-only by default (that is the point — an email per run is noise); `--to ab76@sanger.ac.uk`
adds the email too. It reuses `notify.py` for the credential handling and the `uvx` fallback
below, which is why `notify.py`'s `--to` is optional.

Three things it is careful about, worth not undoing:

- **It always exits 0.** It runs inside the LSF driver job, after the pipeline, and `exit $status`
  in the bsub re-raises Nextflow's real status. A missing `SLACK_BOT_TOKEN`, a node without `lfs`,
  or a `notify.py` traceback must not turn a good run into a failed job — each degrades to
  "unknown" in the message instead.
- **It says the exit status is not the verdict.** `errorStrategy 'ignore'` keeps task failures out
  of the pipeline's exit status (`CLAUDE.md`), so exit 0 is not "nothing failed" — the message
  carries that caveat so a green tick is never read as a clean run.
- **Run name, duration and Nextflow's own verdict come from `.nextflow/history`**, not from the
  log: a killed run leaves `-` in the duration and status columns, which the message reports as
  "killed?" rather than as a clean finish. `lfs quota -q` numbers are KiB and its used value can
  carry a `*` (over soft quota, in grace); the soft limit is preferred over the hard one when set,
  and on scratch124 soft is `0` = unset, so the 142T hard limit is the real ceiling.

## §1 Save the LSF log now, don't wait for the full archive

```bash
.claude/skills/reprocess-debug/bin/save_lsf_log.sh --batch 7
```

`archive_run.sh` (§8) already copies `logs/lsf/reprocessOutput<jobid>.log` and its error twin
into `/nfs/cellgeni/reprocessing-runs/batch<N>/<run>/` — but only as part of archiving the whole
run, which happens "when you are done" with a batch and sometimes doesn't happen for a while.
`logs/lsf/` lives on scratch, which is not backed up. `save_lsf_log.sh` does just this one copy,
safe to run immediately and safe to re-run; archiving the run properly later just overwrites
these two files with the same bytes.

## §2 Refresh the progress counter

```bash
.claude/skills/reprocess-debug/bin/track_progress.py --batch 7 --run elegant_lamarr \
  --json /tmp/progress.json
```

Counts, against the master target list at
`/nfs/cellgeni/projects/reprocessing/irods_datasets/Not_done_hs_10x.sample_table.tsv` (override
with `--target` if a newer snapshot exists — this one is a dated, static export, not something
that refreshes itself), how many of its ~44k samples are done: present on iRODS under
`/archive/cellgeni/datasets/<dataset>/<sample_id>` (any pipeline's upload, not just this repo's —
the target list is project-wide, not scoped to this repo), or completed locally in this repo's
own `results/*/starsolo/<dataset>/<sample_id>/` but not yet uploaded.

**Read this before trusting the number**, every time: the join is by `sample_id` string
equality alone. The target list mixes GEO (`GSM…`), DDBJ (`DRS…`), and other archives' sample
ids in the same column, and this repo's pipeline (`CLAUDE.md`: GEO, SRA, ENA, ArrayExpress) may
never touch a large fraction of the non-GEO rows at all — a low or flat percentage does not mean
this repo is behind, it can just mean most of the list is out of this pipeline's scope. Say that
in the notification rather than letting the number read as this-repo's completion rate.

**The iRODS side reads from a cache**, `/nfs/cellgeni/reprocessing-runs/.cache/irods-datasets-collections.txt`,
refreshed automatically if it is more than 24h old (`--refresh-irods` forces it, `--no-irods`
skips the check entirely if `iquest` isn't available). Regenerating it is one `iquest` call over
~360k collections and takes a few minutes — expected, not a hang. It needs, once per shell,
things that are environment setup rather than this script's job:

```bash
module load cellgen/irods
mkdir -p "$TMPDIR/$USER/singularity"   # or wherever $SINGULARITY_TMPDIR points —
                                        # the icommands wrapper is a Singularity image and
                                        # fails with "could not create temporary directory"
                                        # without this on a fresh shell
```

**A failed `iquest` is recorded as `0`, not as a failure, and the row is appended anyway.** On
2026-09-18 the iRODS session was unauthenticated — `CAT_INVALID_USER`, `iquest` exit 3 — so
`done_irods` came back `0` and batch 9's row reads `1453 / 44494 (3.27%)` with a delta of
**−7833** against batch 8's 20.87%. Nothing in `progress-counter.tsv` marks it as a bad
measurement; the schema has no note column, so the row will read as a real collapse to whoever
opens the file next. The script does print `iRODS : iquest failed (exit 3)` above the numbers —
**read that line before quoting anything**, and if it is there, report the local count only and
cite the last good row instead of the new percentage. `iinit` first (it needs a password, so it
is the user's to run, not something to go looking for) and re-run, or pass `--no-irods` to skip
the iRODS side deliberately rather than silently.

Every run **appends one row** to `/nfs/cellgeni/reprocessing-runs/progress-counter.tsv` — never
overwritten, so the history of the project's progress is the file itself. The script reads the
previous row before appending and reports the delta; that is what makes "how much progress since
last time" a real number rather than a guess.

## §3 Send the notification

```bash
.claude/skills/reprocess-debug/bin/notify.py \
  --to ab76@sanger.ac.uk \
  --subject "Batch 7 debugged — 9 permanent failures, 19.4% of target done" \
  --body-file /tmp/notify-email.txt \
  --slack --slack-text-file /tmp/notify-slack.txt \
  --slack-env-file /lustre/scratch124/cellgen/cellgeni/aljes/reprocessing_slack_bot/.env
```

The *content* is yours to compose each time — what matters changes by batch — but the
*structure* is not: use the house format below for both, so a reader can skim six of these in a
row without re-parsing the shape every time. It always carries, in this order: the artifact link
(or "not published this session"), a short summary of what debugging found, and the progress
count with its caveat. Never print the percentage without the caveat next to it.

**Two different bodies, not one reused for both.** Plain email and Slack `mrkdwn` are different
languages — a bare URL in Slack doesn't render as a link, and `*bold*` in an email shows up as
literal asterisks. Write `--body-file` and `--slack-text-file` as separate files from the
templates below rather than passing the same text to both.

### Email — plain text, section headers, no column alignment

Indentation and a divider line under the title are enough structure; don't try to align numbers
into columns with spaces — most mail clients render plain text in a variable-width font, so
aligned-looking source turns ragged on screen. `label: value` one per line is font-proof. A
section-header emoji is fine (plain text carries Unicode emoji natively) — one per heading, not
one per line; the per-line style below is Slack's, not this one's.

```
🧬 BATCH <N> — <run> — debugged
================================================================

🔗 Post-mortem: <url, or "not published this session">

📋 SUMMARY
- <one line: failed/permanent counts, zero-datasets-lost or not> 🎉 (if zero lost)
- <one line per headline cause, 2-4 total>

🔍 FINDINGS
- 🐛 <n> of <permanent total> — <dataset/accession> — <one-line cause>
  🔧 fix drafted, not applied: issues/<file>.md (if there is one)
- ✅ <n> of <permanent total> — <dataset/accession> — <one-line cause of a correct rejection>

📊 PROJECT PROGRESS  (vs <target file basename>, <total> samples)
- ✅ Done:      <done_total> / <total>  (<pct>%)
    ☁️  on iRODS:                <n>
    💾  local, not yet uploaded: <n>
- ⏳ Remaining: <remaining>
- 📈 Change since last measurement: <+/-n, or "n/a — first measurement">

⚠️ Caveat: <the one-liner from §2 about scope, whenever the percentage is quoted>

--
sent by the reprocess-debug skill's notify step 🤖
```

Pick 🐛 for a real defect and ✅ for a correct rejection consistently — that pairing is what lets
a reader tell the two apart from the emoji alone, which is the point of using them here rather
than decoration for its own sake.

### Slack — `mrkdwn`, not the same text as the email

Slack renders `*bold*`, `_italic_`, `` `code` ``, and `<url|text>` links; it does not render
literal asterisks as anything but asterisks in a client that isn't Slack, which is the other
reason not to reuse this for the email body.

```
🧬 *Batch <N> debugged* — `<run>`
🔗 <url|Post-mortem>

📋 *Summary*
⚠️ <one line, failed/permanent counts> (🎉 if zero datasets lost)
🐛 <n>/<total> → `<dataset>`: <one-line cause> (🔧 fix drafted, not applied)
✅ <n>/<total> → `<accession>`: correct rejection — <one-line reason>

📊 *Progress* (vs `<target file basename>`)
✅ Done: *<done_total> / <total>* (*<pct>%*) — ☁️ <n> on iRODS + 💾 <n> local-only
⏳ Remaining: <remaining>
📈 Δ since last: <+/-n, or "n/a (first measurement)">
⚠️ _Caveat: <the same one-liner, in italics>_
```

Same pairing rule as the email: 🐛 marks an actual defect, ✅ marks a correct rejection, so the
two never blur together at a glance. `reprocessing_slack_bot/app.py`'s own style
(`✅❌:among_us_hammer::bird_run:`) is that bot's own custom-emoji workspace shorthand for a
different message shape (pass/fail/processed/left counts) — borrow the *idea* of a symbol per
state, not those specific glyphs, since `:among_us_hammer:` and `:bird_run:` are workspace emoji
that may not exist wherever this posts.

**Email needs no setup** — `smtplib.SMTP("localhost")` uses the farm's local relay, same as
`sample-tracking`'s `--email`. It will simply work.

**Slack needs `SLACK_BOT_TOKEN` and `CHANNEL_ID`**, read from the environment or from
`--slack-env-file` (a plain `KEY=VALUE` file, e.g. `reprocessing_slack_bot/.env` — a *different*
project's file, reused here only for its credentials, not imported as a dependency). If neither
resolves to a real value, `notify.py` reports "skipped" on stderr and still sends the email —
**never fabricate or go looking for a token elsewhere**; if Slack matters, tell the user which
two variables are needed and where to put them, and let them decide.

**Slack's library comes from `uvx`, and `pip install --user` cannot substitute for it.** The farm
python has user site-packages **disabled** (`site.ENABLE_USER_SITE` is `False`), so
`pip install --user slack_sdk` reports success, writes
`~/.local/lib/python3.10/site-packages/slack_sdk`, and nothing ever imports it — which is exactly
how batch 8's notification went out with `slack: slack_sdk not installed in this environment —
skipped` on 2026-09-16. `notify.py` now tries the in-process import first and, on `ImportError`,
re-runs the one API call inside an ephemeral environment:

```bash
uvx --with slack_sdk python -c '<the post>'      # or `uv run --no-project --with slack_sdk …`
```

`uv` and `uvx` are already on `PATH` from `/software/cellgen/cellgeni/uv/` (module
`cellgen/uv/0.11.8` is the default if they are not). The success line says which route was used —
`slack: posted to channel <id> (slack_sdk via uvx)` — so "posted" without that suffix means the
module was importable directly. The token and channel reach the child through its **environment,
not argv**: `/proc/<pid>/cmdline` is world-readable on a shared login node and
`/proc/<pid>/environ` is not. If neither `uvx` nor `uv` is on `PATH`, Slack skips with that said
plainly and the email still sends.

The first `uvx` call in a session downloads `slack_sdk` (~320 KB) and takes a second or two;
later calls are cached.

**Tell the user, in the chat reply, whether Slack actually sent** — `notify.py`'s stderr says so
plainly (`slack: posted to channel …` vs `slack: … — skipped`); don't imply it went out unless
that line confirms it did.
