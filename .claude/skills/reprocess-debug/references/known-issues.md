# Baseline distributions and open issues

Distilled from `docs/agent_debug.md` §7 and §13, plus what the September 2026 collections added.

## §baseline — the August 2026 run (1814 failures)

A reference distribution. If a new run's profile differs sharply from this, that difference is
itself the finding.

Stage: `RENAME10XRUN` 1738 · `RENAME10XSAMPLE` 39 · `STARSOLO10X` 23 · `WGET10X` 14.
Totals: 1800 classified + 14 with no log block = 1814.

| Cause | Runs | Verdict |
|---|---|---|
| `chem-no-whitelist-hit` | 1014 | 880 genuinely not 10x (GSE109816); the rest mixed |
| `bc-variable-length` | 401 | pipeline too strict — chemistry called at 96%+ |
| `no-single-gex-r2` | 159 | correct: feature-barcode (ADT/HTO) libraries, R2 ≈ 25 bp |
| `single-fastq-from-sra` | 89 | `SRA2FASTQ` emitted one mate |
| `read-counts-differ` | 48 | `SRA2FASTQ` partial dump |
| `r1-len-differs-runs` | 39 | missing chemistry JSON wiring |
| `chem-no-whitelist-rand` | 13 | STARsolo re-detect disagreeing with run-level inference |
| `chem-ambiguous-layout` | 13 | too strict |
| `truncated-gzip` | 9 | download integrity |
| `star-sj-buffer` | 6 | `--limitOutSJcollapsed`; 5 at exit 104, 1 at 139 |
| `corrupt-fastq` | 5 | malformed record mid-file |
| `UNCLASSIFIED` | 4 | 3 = LSF `TERM_MEMLIMIT` at exit 130; 1 = a `starsolo_10x_auto.sh` check |
| no log block | 14 | `WGET10X`: 13 × HTTP 403 (transient), 1 × empty URL |

Roughly 1040 correct rejections, 774 avoidable.

## §later — what the later batches looked like

Failure volume dropped by two orders of magnitude once the guards landed. The remaining
failures are a different population, so do not expect the August shape.

| Batch | Run | Failed tasks | Permanent | Notes |
|---|---|---|---|---|
| 3 | `big_keller` | 111 | 92 | 41 UNCLASSIFIED — that manifest predates the LSF-kill line, so exit 130/139 carry no evidence |
| 5 | `tender_brattain` | 48 | 48 | 32 `sra-zero-length-read` in one dataset (GSE200629 = 24 of 48) |
| 6 | `spontaneous_ampere` | 16 | 7 | 9 self-healed on retry; 3 × `chem-no-whitelist-rand`, 2 × `star-segfault`, 1 × `meta-no-run-id`, 1 × `sra-too-few-fastq` |

Batch 6 is the useful worked example: `.nextflow/history` records it `OK`, and it still had 7
permanent failures. `errorStrategy 'ignore'` means a green run tells you nothing.

## §open — found but not fixed

* **`errorStrategy` precedence bug.** `&&` binds tighter than `||`, so every task retries once
  including 128 GB STAR jobs on deterministic errors. `attempt` is not diagnostic.
* **No failure manifest inside the pipeline.** `bin/collect_run_logs.sh` reconstructs one after
  the fact, but a `workflow.onComplete` hook writing run/sample/dataset/process/exit/first-error
  would remove the whole collection step. A previous attempt at this was reverted (commit
  `028aaf3`) for causing silent failures.
* **Chemistry is still detected twice.** `chemistry.json` reaches `RENAME10XSAMPLE`, but
  `starsolo_10x_auto.sh` re-detects from scratch and can disagree with run-level inference. This
  is `chem-no-whitelist-rand`: 13 late failures in August 2026, still 3 in batch 6.
* **Empty `chemistry.json` reaches `RENAME10XSAMPLE`.** Batch 5's two exit-1
  `RENAME10XSAMPLE` failures are a bare `JSONDecodeError` out of `load_run_metadata` — the file
  exists and is empty. The classifier calls this `chem-json-unreadable`; there is no guard.
* **`FETCH10XMETA` can fail per-sample with no run ID.** `No experiment or run ID found for
  GSM… in GSE….sra.tsv` (batch 6, GSE206528). Not in the August taxonomy at all.
* **No MD5 verification.** ENA publishes `fastq_md5`, but plumbing it through means extending
  `curl_ena_metadata.sh`, the awk join in `parse_metadata.sh`, `collect_metadata.sh` and the
  searchlist schema — and it only covers the ENAFQ route. `gzip -t` catches truncation but not a
  valid gzip holding malformed FASTQ.
* **No pre-download 10x screen.** GSE109816 (880 runs, not 10x) was downloaded in full before
  being rejected at `RENAME10XRUN`. A per-dataset whitelist sample right after `FETCH10XMETA`
  would stop that.
* **Feature-barcode libraries are not filtered on metadata.** ADT/HTO/hashing libraries are
  usually labelled as such in GEO and could be skipped before download.
* **`REPROCESS10X_BAM2FASTQ` stub is broken** — it touches `<id>.fastq.gz` but declares
  `path("fastqs/*")`, blocking stub-testing of any BAM-origin path.
* **`modules/local/reprocess10x/sra2fastq/` is dead code** — the more careful implementation,
  not wired up. Wire it up or delete it.
* **`scripts/run_reprocess.bsub` does not pass `-name`.** Every run gets a generated name and
  the batch → run link has to be recovered from `.nextflow/history`. One flag would fix it.
* **Module drift.** The bsub loads `cellgen/nextflow/26.04.1`, CLAUDE.md and README say
  `26.04.6`, the manifest requires `>=26.04.1`, and the PATH default is 25.04.4 — which reads a
  26.x `.nextflow/cache` as an empty run and rejects the `accelerator` trace field.
