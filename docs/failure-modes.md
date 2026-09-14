# Failure modes

What goes wrong when this pipeline runs at scale, what each failure looks like, and whether the
cause is the pipeline or the data.

This is the orientation page. The debugging *procedure* — collecting a run's logs, classifying
them, verifying a claim against the work directories — lives in the `reprocess-debug` skill
under [`.claude/skills/reprocess-debug/`](../.claude/skills/reprocess-debug/), which is the
right place to start if you have a specific run to triage.

## Two things that make failures hard to read

**A run that finished `OK` tells you nothing.** Every process sets `errorStrategy 'ignore'`, so
failed tasks never reach the pipeline's exit status and there is no failure manifest. Batch 6 is
recorded `OK` in `.nextflow/history` and had 7 permanent failures. Always count the failures;
never infer from the status.

**`attempt` is not diagnostic either.** The `errorStrategy` in
[`nextflow.config`](../nextflow.config) gives *every* task one retry whatever its exit code, and
only the codes in the retriable set any more than that. That is intentional — see
[Which failures are worth retrying](#which-failures-are-worth-retrying). The effect on triage is
that in August 2026, 1795 of 1814 failures sat at `attempt = 2`: a property of the config, not a
signal about any of them. A task that failed twice is not more interesting than one that failed
once; what matters is whether any attempt later succeeded.

Because of the retry, always separate **permanent** failures from ones that self-healed. Batch 6
had 16 failed tasks, 9 of which recovered on retry — the real number is 7, and every percentage
should be over that.

## Exit codes

| Process | Exit | Meaning |
|---|---|---|
| `RENAME10XRUN` | 2 | `TenxRunError` from the inference script, **or** an argparse usage error — argparse dumps the whole option list, so read the `error:` line at the end |
| `RENAME10XSAMPLE` | 2 | `RenameError` from the rename script |
| `RENAME10XSAMPLE` | 1 | Uncaught Python exception — an empty `chemistry.json` surfaces here as a bare `JSONDecodeError` |
| `WGET10X` | 8 | An HTTP error response (403/404). `wget` treats these as fatal and ignores `--tries` unless `--retry-on-http-error` is set |
| `WGET10X` | 1 | `wget` usage error — in practice an empty URL field in the metadata |
| `SRA2FASTQ` | 1 | A post-dump validation guard tripped: wrong number of FASTQs, mates at different depths, a read discarded as zero length |
| `STARSOLO10X` | 1 | A check inside `starsolo_10x_auto.sh`, e.g. `No whitelist matched 200,000 random barcodes` |
| `STARSOLO10X` | 104, 139 | STAR's own fatal error. 104 has been the splice-junction buffer (`--limitOutSJcollapsed`); 139 is SIGSEGV |
| any | 130 | 128+2, but in practice an LSF `TERM_MEMLIMIT` kill, not a signal from the tool. Confirm against the `TERM_*` line in the raw log before treating it as a crash |
| any | 140, 146 | Other LSF kills — `TERM_RUNLIMIT`, queue eviction |

Exit 130 is the one that misleads. It looks like an interrupt and is almost always the job
running out of memory.

## Which failures are worth retrying

Every task gets one retry whatever its exit code; the codes below get a third attempt. The
`task.attempt == 1` clause in [`nextflow.config`](../nextflow.config) is deliberate and
load-bearing — **do not remove it.** An exit code does not tell you whether a failure was
transient, so everything gets a second chance, and the cost of rerunning a doomed task is
accepted in exchange for not discarding a recoverable one. A task is also only ignored after
failing *twice*, which is what keeps the first attempt's work directory and stderr around for
triage.

Retried up to three times: **130-145** (the signal range — 130 is in practice an LSF
`TERM_MEMLIMIT` kill, 139 a SIGSEGV), **104** (STAR's own fatal error), **175** (LSF requeue).

What has actually been observed, from batch 6 — 16 failed tasks, 9 of which recovered:

| Exit | First error | Attempts failed | Recovered | Read as |
|---|---|---|---|---|
| 130 | `[lsf] TERM_MEMLIMIT` | 1 | 6 of 6 | transient — the node was oversubscribed |
| 104 | `FATAL ERROR in reads input: quality string length` | 1 | 1 of 1 | transient — read off a FASTQ still being staged |
| 1 | `R1 length (-2147483647) leaves no room for a UMI` | 1 | 1 of 1 | transient — a nonsense length from the same cause |
| 139 | Segfault in STAR under `--genomeLoad LoadAndRemove` | 3 | 0 of 2 | deterministic here, though the shared-memory genome is a genuine race |
| 1 | `No whitelist matched 200,000 random barcodes` | 2 | 0 of 3 | deterministic — the data is not what we think it is |

**These are counts from one batch, not a law.** Row 3 is the one that justifies the design: exit
1 is not in the retriable set and never should be, yet that sample was recovered purely by the
blanket retry. Rows 4-5 are what it costs.

When you identify a new transient class, add its exit code to the set in `nextflow.config` *and*
a row to this table in the same change — otherwise the next person rediscovers it from scratch.

## Where failures come from

A reference distribution, from the August 2026 run of 1814 failures. A new run whose shape
differs sharply from this is itself worth investigating.

Stage: `RENAME10XRUN` 1738 · `RENAME10XSAMPLE` 39 · `STARSOLO10X` 23 · `WGET10X` 14.

| Cause | Runs | Whose fault |
|---|---|---|
| No whitelist hit at all | 1014 | Mostly the input list — 880 were GSE109816, which is not 10x |
| Barcode length varies within a run | 401 | Pipeline: too strict, chemistry was called at 96%+ confidence |
| No single GEX R2 | 159 | Correct rejection — feature-barcode (ADT/HTO) libraries, R2 around 25 bp |
| `SRA2FASTQ` emitted one mate | 89 | Archive |
| Mates dumped at different depths | 48 | Archive |
| R1 length differs between runs of a sample | 39 | Pipeline: missing chemistry wiring |
| STARsolo re-detect disagrees with run-level inference | 13 | Pipeline: chemistry detected twice |
| Ambiguous layout | 13 | Pipeline: too strict |
| Truncated gzip | 9 | Download integrity |
| STAR splice-junction buffer | 6 | Pipeline: needs `--limitOutSJcollapsed` |
| Corrupt FASTQ mid-file | 5 | Archive |

Roughly 1040 correct rejections against 774 avoidable ones. Failures also concentrate hard —
the top five datasets were 82% of the total and one alone was 48% — so **aggregate by dataset
before drawing any conclusion.** A flat list of failed accessions hides the concentration, and a
fix aimed at the flat list aims at the wrong thing.

Volume dropped by two orders of magnitude once the guards landed; later batches are a different
population and should not be compared to these proportions directly.

## Open problems

Known, reproduced, and not yet fixed.

**Data is downloaded before it is screened.** Nothing checks `library_strategy` between metadata
fetch and download, so correctly labelled ATAC-seq, bisulfite-seq and bulk RNA-seq are fetched
in full and rejected afterwards — 596.8 GB in one run. The metadata needed is already on disk
when the decision is made. This is the single largest saving available and is listed first for
that reason. See [Archive pathologies](archive-pathologies.md#correctly-labelled-rejected-late).

**Chemistry is detected twice and the two can disagree.** `RENAME10XRUN` writes a per-run
`chemistry.json`, but `starsolo_10x_auto.sh` re-detects from scratch and can reach a different
answer. Separately, Cell Ranger auto-detects ARC-v1 for Multiome GEX libraries and then refuses
to proceed without being told explicitly — costing the GEX halves of two datasets that the
pipeline had already correctly identified upstream.

**An empty `chemistry.json` reaches `RENAME10XSAMPLE`** and surfaces as a bare `JSONDecodeError`.
The file exists and is zero-length; there is no guard.

**Run-to-sample mapping is keyed on BioSample.** Where a submitter registers two libraries
against one BioSample, the last writer wins and two libraries merge silently into one sample.
Keying on the experiment (SRX) instead, and failing loudly when a run resolves to more than one
sample, would close it. See the GSE247111 row in
[Archive pathologies](archive-pathologies.md#defects).

**No checksum verification.** ENA publishes `fastq_md5` and it is not used. `gzip -t` catches
truncation but not a valid gzip containing a malformed FASTQ.

**No failure manifest.** Failures have to be reconstructed after the fact from the Nextflow
trace and the work directories. An earlier attempt at a `workflow.onComplete` hook was reverted
for causing silent failures.

**Runs are not named.** The launch script does not pass `-name`, so every run gets a generated
`adjective_surname` and the batch-to-run link has to be recovered from `.nextflow/history`. One
flag would fix it.
