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

**`attempt` is not diagnostic either.** The `errorStrategy` expression in
[`nextflow.config`](../nextflow.config) has an operator-precedence bug — `&&` binds tighter than
`||`, so `... && task.attempt < 3 || task.attempt == 1` retries *every* task exactly once
whatever the exit code, including 128 GB STAR jobs failing deterministically. In August 2026,
1795 of 1814 failures sat at `attempt = 2`, which was a property of the config and not a signal
about any of them. A task that failed twice is not more interesting than one that failed once;
what matters is whether any attempt later succeeded.

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
