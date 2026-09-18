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
| 7 | `elegant_lamarr` | 10 | 9 | 1 self-healed (`lsf-memlimit`); 8 × `chem-ambiguous-layout` (one dataset, GSE212038 — truncated-UMI allowlist gap, see `§open`), 1 × `chem-no-whitelist-hit` (genuinely near-empty source FASTQ, confirmed against the raw download). First batch with zero datasets lost. |
| 8 | `cheeky_cuvier` | 53 | 11 | 42 self-healed (35 × `wget-exhausted`, 6 × `truncated-gzip`, 1 × `lsf-memlimit`); **10 of the 11 permanent are one pathology in one dataset pair** — GSE212964/GSE212965 split-mate runs, 5 samples, see `§open` — plus 1 near-empty SRA object (GSE213370). 6 of 787 unique samples lost. Also the first batch where duplicate dataset rows were quantified: 200 of 981 alignment tasks redundant (20.4%). |
| 9 | `gigantic_hypatia` | 108 | 102 | 6 self-healed (3 × `truncated-gzip`). **The failure list is the small half of this batch.** 68 of the 102 are correct rejections (54 ADT/CRISPR runs in GSE215253 alone); the two real defects are 9 samples lost to a STAR segfault in GSE214532 (`§open`) and 2 to a `umi_len` disagreement in GSE215362. Separately, and invisible in the failure list: **38 aligned samples are not GEX libraries at all** and were published as near-empty matrices, and 133 of 917 alignments (14.5%, 2131 CPU-h) were duplicates. 55 of 839 samples lost. |

Batch 6 is the useful worked example: `.nextflow/history` records it `OK`, and it still had 7
permanent failures. `errorStrategy 'ignore'` means a green run tells you nothing.

## §open — found but not fixed

* ~~**`errorStrategy` precedence bug.**~~ **Closed — not a bug.** Retrying every task once
  whatever its exit code is intentional, and the `task.attempt == 1` clause is load-bearing:
  batch 6 lost no sample to a transient exit-1 failure because of it. See `classify.md §exits`
  for the evidence table. Two things remain true and are not defects: `attempt` is not
  diagnostic, and a deterministic failure costs one extra attempt — the accepted price of not
  discarding transient ones.
* **No failure manifest inside the pipeline.** `bin/collect_run_logs.sh` reconstructs one after
  the fact, but a `workflow.onComplete` hook writing run/sample/dataset/process/exit/first-error
  would remove the whole collection step. A previous attempt at this was reverted (commit
  `028aaf3`) for causing silent failures.
* **`--allow-gex-truncated-umi-lengths` default (`10`) is missing `11`.** In
  `infer_10x_run_recommended.py`, this flag is the guard's own named mechanism for a v3 GEX read
  that is short on UMI bases, but its default only accepts an observed UMI of 10. Batch 7's
  GSE212038 (8 runs, `SRR21201139/156/157/158/159/164/165/166`) has R1 uniformly 27 bp — 16 bp CB
  + 11 bp UMI — which falls through the allowlist and gets the generic "no supported 10x run
  layout matched" rejection instead of the truncated-UMI acceptance path. Confirmed by re-running
  the script against the persisted work dir with `--allow-gex-truncated-umi-lengths 10,11`: it
  then accepts the run with a truncated-UMI warning, no other change. Fix is one value, in one
  place — `ext.args` for `RENAME10XRUN` in `configs/rename10xrun.config`. Not yet applied; see
  the [Batch 7 Post-Mortem](../../../../docs/post-mortems.md) for the full evidence.
* **Split-mate runs are rejected per run, upstream of the merge that would fix them.** When a
  submitter registers R1 and R2 of one 10x library as two SRA run accessions, `RENAME10XRUN`
  infers layout per run and rejects both halves — the barcode half for having no biological
  read, the cDNA half for having no whitelist hit — although `links.tsv` already names the same
  GSM on both rows and `RENAME10XSAMPLE` exists downstream to merge them. Cost so far: 5 samples
  in batch 8 (GSE212964/GSE212965, measured) and 1 in batch 5 (GSE202476, log-text signature
  only — its work dirs are gone). The detection signature is a sample whose runs fail with
  *complementary* errors. Fix needs a design decision, not a threshold; drafted in
  `issues/2026-09-16-batch8-split-mate-runs.md`.
* **Datasets with identical sample lists are processed twice.** GEO SuperSeries/SubSeries pairs
  whose `sample_id` field is byte-identical become two dataset units, so every shared sample is
  downloaded, renamed and aligned twice. Batch 8: 200 of 981 alignment tasks redundant
  (20.4% of the stage, 26% more than the 781 distinct samples needed), plus 545 redundant
  `WGET10X` and 56 redundant `SRA2FASTQ`. Batch 9: 21 duplicated pairs of 118 dataset rows, 133
  of 917 alignments (14.5%) and 2131 of 14 844 CPU-hours. Project-wide, 445 of 4602 dataset rows
  (9.7%) are exact duplicates. Not a correctness bug — the outputs are right — but the largest
  single waste measured so far, and it multiplies failures too: GSE214532's segfault burned 54
  attempts rather than 27 because the dataset ran twice. It also means a per-sample tally taken
  off `mapping_qc_stats.tsv` without de-duplicating is over by 17%. Drafted in
  `issues/2026-09-16-batch8-duplicate-dataset-rows.md`.
* **`triage.py` verdicts `no-single-gex-r2` as a correct rejection.** The note asserts an
  ADT/HTO library, but the split-mate barcode half emits the same error text, and in batch 8 all
  5 rows were the latter — so the triage output read 54.5% correct-rejection for a set that was
  really 1 correct rejection and 10 tasks of one defect. Batch 5's 7 rows were misverdicted the
  same way (4 of them the GSE202476 sample). One-line fix, dry-run confirmed, drafted in
  `issues/2026-09-16-batch8-triage-verdict-split-mate.md`.
* **`collect_run_logs.sh` writes a ragged trace when `nextflow log` warns on stdout.** Nextflow
  puts its FTP-proxy warnings on stdout, which is redirected straight into `runlogs<N>.tsv`; one
  stray line and duckdb's dialect sniffer refuses the whole file, with an error that mentions
  only CSV dialects. It bites from an agent session, whose sandbox sets `FTP_PROXY` to a
  `socks5h://` URL that Nextflow will not parse — i.e. exactly how this skill is normally run.
  All surviving traces (batches 6-8) are clean at 36 fields. Drafted, with diff, in
  `issues/2026-09-16-batch8-collect-trace-ragged.md`.
* **Nothing gates a sample on whether it mapped.** Batch 9 published 38 near-empty matrices —
  HTO/hashing, ADT/CSP, CRISPR-enrichment and "custom" feature libraries aligned as GEX. Every
  task completed, so none of it is in `failures9.tsv`; the pathology is only visible in
  `mapping_qc_stats.tsv`, and only if you compute it. Neither metadata nor geometry can filter
  these: all 38 are `library_strategy = RNA-Seq` / `transcriptomic single cell`, and GSE216914's
  CSP and VDJ libraries share an identical `10/10/26/90` layout while mapping at 0.9% and 93.9%.
  The wrapper's own strand test already measures the discriminating number and warns without
  acting — a gate at `max(fwd,rev) < 5%` flags 31 samples with zero false alarms, and post hoc
  `exon_u < 0.05` separates the set completely (flagged max 0.0398, kept min 0.0597). The gate
  belongs in the container's `starsolo_10x_auto.sh`; a repo-side `Flag` column for
  `mapping_qc_stats.tsv` is drafted, with diff, in
  `issues/2026-09-18-batch9-no-mapping-rate-gate.md`. **Do not chase this with a read-length
  threshold** — 42 bp of feature library is junk and 55 bp of GEX maps at 94%, measured in the
  same batch.
* **The STARsolo strand test decides paired-end mode on no margin, and the paired-end STAR call
  segfaults.** `starsolo_10x_auto.sh` sets `PAIRED=True` from `R1LEN > 50` alone, then keeps it
  unless the strand test returns `Forward`; the test concludes `Reverse` on `PCTREV > PCTFWD`,
  strictly greater, with no floor. Samples that stay paired-end run
  `--soloBarcodeMate 1 --clip5pNbases 39 0` and STAR 2.7.10a_alpha_220818 segfaults immediately
  after thread creation. Batch 9's GSE214532: 12 samples, identical 150+150 geometry, 9 lost at
  exit 139 on 6 attempts each (54 attempts, 18 CPU-h), 3 completed at 94-95% mapping — one of
  them on a 44/44 tie. All 12 tripped the script's own "Low percentage of reads mapping to
  GeneFull" warning, so the test was inconclusive for every one of them and was still allowed to
  decide. The file is inside `quay.io/cellgeni/starsolo:v4.3` and there is no flag to override it
  from here; drafted in `issues/2026-09-18-batch9-starsolo-strand-test-segfault.md`. Batch 6's 2
  `star-segfault` rows and batch 3's 1 may be the same branch — unverifiable, their work dirs are
  gone.
* **`umi_len` is held to exact equality across a sample's runs, and it is an observed quantity.**
  A run sequenced with a longer R1 reports 12 observed UMI bases where its siblings report 10, and
  `require_metadata_compatibility` rejects the sample. Batch 9 lost GSM6634355 and GSM6634358
  (GSE215362, 9 runs each, 8-vs-1) this way — same chemistry id, same CB length, whitelist match
  93.7% and 95.8%. The run-level inference already tolerates the short case, warning "treating as
  truncated-UMI layout", and then the sample-level guard rejects exactly what it accepted. Fix is
  to take the minimum; drafted in `issues/2026-09-18-batch9-umi-len-mismatch-rejects-sample.md`
  (no diff — registry module with a `.module-info` checksum).
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
