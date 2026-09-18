# Drafting issues for confirmed findings

Every debugging session that turns up a confirmed, actionable defect drafts an issue for it —
not just runs that also get a post-mortem. A post-mortem is for a person; an issue is for the
tracker, and the two are written from the same evidence but are not substitutes for each other.

**Never call `gh issue create` and never apply the suggested patch to the working tree.** Both
files in this section are drafts for the user to review — see `SKILL.md`'s note on this. Filing
the issue and applying the fix are the user's decisions, made after reading what you wrote here.

## §when

Draft one for a finding whose `verdict` (from `triage.py`, or your own `verify.md` measurement
overriding it) is a real defect: `too-strict`, `pipeline-bug`, `investigate` that resolved to a
confirmed cause, or "our tooling" — anything that would otherwise become a new bullet in
`known-issues.md §open`. In fact, if you are about to add such a bullet there, that is the
trigger: draft the issue in the same pass.

Do **not** draft one for a `correct-rejection` verdict (batch 7's `SRR20761283` — a genuinely
near-empty source file — got no issue) or for an `infrastructure`/self-healed row (a
`TERM_MEMLIMIT` that succeeded on retry is not a defect; see SKILL.md rule 1). One issue per
distinct root cause, not per failed task — batch 7's eight `GSE212038` runs are one issue, not
eight.

**Check `./issues/` for an existing draft on the same root cause before writing a new one.**
Batches recur against the same guards; a second batch hitting the same gap updates the existing
file's evidence section rather than forking a duplicate.

## §where

`./issues/` at the repo root, created if it does not exist. Nothing in `.gitignore` catches it or
the `.md`/`.diff` files inside — the extension rules (`*.json`, `*.log`, `*tsv`, …) do not list
either, and `issues` itself is not a named pattern. This is deliberate: these are meant to be
reviewable and eventually committed once the user has approved them, unlike `data/` or
`scripts/`.

**Name every file with today's date first.** Debug sessions stack up in `./issues/` between
reviews — the user works through them in a batch, not one per session — so a batch number alone
is not enough to tell a fresh finding from one that has been sitting for three weeks, and two
sessions on the same day still need `date` to sort before `ls` does anything useful. Get the
date from the environment (`date +%F`), never guess it:
`issues/<YYYY-MM-DD>-batch<N>-<short-slug>.md` and, when there is a concrete fix,
`issues/<YYYY-MM-DD>-batch<N>-<short-slug>.diff` alongside it, same date and slug. Use the date
and batch/run that *found* the issue, not necessarily the only one it affects — if a later batch
hits the same root cause, that updates the existing dated file's evidence section (§when) rather
than taking a new date.

## §the issue file

Plain markdown, written the way it will be pasted into `gh issue create --body-file`. No house
style, no palette — this is not the Artifact report.

```markdown
# <one-line, imperative or descriptive title — this becomes the GH issue title>

**Drafted:** <YYYY-MM-DD, matching the filename — not the batch's run date>
**Found in:** batch <N>, run `<name>` — <link to its post-mortem if one was published>
**Suggested labels:** <e.g. bug, chemistry-inference>
**Suggested fix:** see `issues/<YYYY-MM-DD>-batch<N>-<short-slug>.diff`

## Summary

One paragraph: what breaks, for whom, how often.

## Evidence

The measured facts — accessions, work-dir paths, the specific numbers from `chemistry.json` /
the trace / a direct re-run — exactly as you'd have verified them per `verify.md`. Quote the
tool's own error text once; don't paraphrase it.

## Root cause

File and line. If you found the exact branch or default that causes it, quote it.

## Suggested fix

What changes, and why that's the minimal correct change rather than a broader rewrite. Note any
value picked (a threshold, an allowlist entry) that is a judgment call the user should confirm,
not a fact you measured.

## Regression / testing notes

What to re-run after applying the fix — `verify.md §cases` regression accessions if inference is
touched, `verify.md §regression` if `triage.py`'s rules changed, or a stub run if it's channel
wiring. Say if the fix has been dry-run confirmed already and how.

## Not covered

Anything the fix doesn't address — a related case you didn't check, a broader version of the
same gap you chose not to fix pre-emptively.
```

## §the suggested diff

Generate a real, byte-correct unified diff — do not hand-type one from memory of the file's
content, line numbers drift and a bad diff is worse than no diff. The safe procedure, which
never leaves the working tree modified:

```bash
# 1. Make the fix for real, with Edit — e.g. configs/rename10xrun.config
# 2. Capture it
git diff -- configs/rename10xrun.config > issues/2026-09-15-batch7-truncated-umi-allowlist.diff
# 3. Revert the working tree — the fix lives in the diff file, not in the tree
git checkout -- configs/rename10xrun.config
# 4. Confirm the diff you just wrote actually applies
git apply --check issues/2026-09-15-batch7-truncated-umi-allowlist.diff && echo "applies clean"
```

Multiple files: pass all of them to both the `git diff` and the `git checkout` in steps 2 and 3.
Step 4 is not optional — an issue with a diff that fails `git apply --check` is worse than one
with no diff, because it looks actionable and isn't.

If the fix isn't a small, confident change — it needs a design decision, or touches a shared
module rather than a per-process config — write the issue without a diff and say so in
**Suggested fix** rather than guessing at a bigger patch.

## §telling the user

State plainly, in the chat reply, that both files are drafts pending their review — name the
paths, and say neither `gh issue create` nor `git apply` has been run. This is not a place to be
terse: the whole point of this step is that nothing ships without a human reading it first.
