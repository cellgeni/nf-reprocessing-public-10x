# Run index

Every reprocessing run, what it was, where its evidence is archived, and which post-mortem
covers it. **Check here before hunting** — the artifact URLs in particular are otherwise only
findable by listing every published artifact and guessing from titles.

Archive root: `/nfs/cellgeni/reprocessing-runs/`, one directory per run, each with a
`manifest.txt` describing itself. Post-mortem URLs are `https://claude.ai/code/artifact/<uuid>`
and resolve only for the account that published them.

`.nextflow/history` is the authority on the first five columns; see `data-map.md §history`.

> ### The work directories for batches 1-5 are gone
>
> Deleted, confirmed 2026-09-14. **The archive is the only surviving record for those runs**, and
> no claim about them can be re-verified against a FASTQ any more — `verify.md §walkback` and
> SKILL.md rule 3 ("measure, do not infer") apply only to batch 6 onward.
>
> Batch 4 shows what that costs. It was collected for the first time on 2026-09-14: the Nextflow
> trace survived, so we know it had **47 failed tasks, 44 permanent** — but every one of the 47
> carries `first_error = (no .command.log found)`, because the stderr only ever existed in the
> work directory. The counts are recoverable from the trace; the reasons are not.
>
> **Collect and archive a run while its work dirs still exist.**

## Batch runs

| Batch | Run | Started | NF status | Session | Failed | Permanent | Archive | Post-mortem |
|---|---|---|---|---|---|---|---|---|
| 1 | `peaceful_feynman` | 2026-08-19 13:27 | OK | `9fd63dba` | — | — | `batch1/peaceful_feynman` | — |
| 2 | `distraught_mercator` | 2026-08-21 15:13 | OK | `9fd63dba` | 294 | unknown † | `batch2/distraught_mercator` | — |
| 3 | `big_keller` | 2026-08-25 18:24 | OK | `9fd63dba` | 73 in log † | 92 of 111 ‡ | `batch3/big_keller` | [Batch 3 Failure Triage](https://claude.ai/code/artifact/ef2992f4-b281-4942-8f17-8ceb4c90fbb5) |
| 4 | `nauseous_spence` | 2026-09-01 14:55 | OK | `4645283e` | 47 | 44 | `batch4/nauseous_spence` | — (no stderr; see box above) |
| 5 | `tender_brattain` | 2026-09-08 12:43 | OK | `f97ddf5e` | 48 | 48 ‡ | `batch5/tender_brattain` | [Batch 5 Post-Mortem](https://claude.ai/code/artifact/9d667def-0edb-4278-9a6d-50bc49c92931) |
| 6 | `spontaneous_ampere` | 2026-09-08 23:25 | OK | `34653a63` | 16 | **7** | `batch6/spontaneous_ampere` | — |

**Every one of these is recorded `OK`.** `errorStrategy 'ignore'` keeps failed tasks out of the
exit status, so the status column is not a verdict — batch 6 is `OK` with 7 permanent failures.

† **`Permanent` is only meaningful where a `failures<N>.tsv` manifest exists, and only batch 6
has one.** For the others `triage.py` was run over the `failed<N>.log` / `failedjobs<N>.tsv`
pair, which cannot know which tasks later succeeded on retry, so its `permanent` field equals
its `total`. Do not quote those totals as permanent-failure counts.

‡ From the analysis recorded in `known-issues.md §later`, which reconciled recovery properly.
Batch 3's 111 failed tasks / 92 permanent comes from there; the 73 in the archived
`triage3.json` is the number of blocks in `failed3.log`, which is smaller. Prefer the
`known-issues.md` figures for batch 3 and 5.

### Runs superseded by a resume

A `-resume` keeps the session uuid and takes a new run name. Triaging the wrong one of a pair
shows failures a later resume fixed, or hides them.

| Batch | Superseded run | Started | Replaced by |
|---|---|---|---|
| 4 | `fervent_raman` (ERR) | 2026-08-26 10:41 | `nauseous_spence` — different session, so a fresh attempt rather than a resume |
| 5 | `awesome_neumann` (ERR) | 2026-09-07 11:18 | `tender_brattain`, same session `f97ddf5e` |

## Not batch runs

| What | Archive | Post-mortem |
|---|---|---|
| **August 2026 baseline** — 1814 failed tasks, the reference distribution quoted throughout the docs. Collected by hand on 2026-08-17 and **not attributable to a named run**: it predates `peaceful_feynman` by two days. See its `manifest.txt`. | `legacy/august-2026-baseline` | [Reprocessing Run Post-Mortem](https://claude.ai/code/artifact/e4947587-de6b-4793-bd66-23da91af549a) |
| **REQ-74217** — 21 datasets, mixed public and internal. Source of the GSE247111 shared-BioSample finding and the ARC-v1 chemistry bug. | not archived here | [REQ-74217 Post-Mortem](https://claude.ai/code/artifact/bcb35389-362e-4e89-9309-597e073aa66b) |
| **Run 337224** (LSF job id; Nextflow run `evil_volhard`) — the REQ-74217 follow-up that found the GSE247111 library merge. | not archived here | [Run 337224 Post-Mortem](https://claude.ai/code/artifact/7d89d7e0-705a-4c5b-bdb3-621615a9504c) |
| **Run 673887** (LSF job id) — GSE188429, GSE189141, GSE183068, GSE188823. | not archived here | [Run 673887 Post-Mortem](https://claude.ai/code/artifact/547e62d8-45bb-4aa6-ad6a-0a95cd9d5f36) |

The last three ran in a **different working directory**: neither `evil_volhard` nor LSF jobs
337224 and 673887 appear in this repo's `.nextflow/history` or `logs/lsf/`, so their traces and
work dirs are not reachable from here. The post-mortems are the only surviving record.

## Adding a run

```bash
.claude/skills/reprocess-debug/bin/collect_run_logs.sh --batch 7
.claude/skills/reprocess-debug/bin/archive_run.sh      --batch 7
```

Then add the row. When a post-mortem is published, put its URL in the last column — that is the
whole point of this file, and it is the step most easily forgotten because the artifact is
already safely stored and feels done.
