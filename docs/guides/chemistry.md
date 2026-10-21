---
title: Chemistry inference
parent: Guides
nav_order: 4
---

# Chemistry inference
{: .no_toc }

1. TOC
{:toc}

## Why the pipeline infers it

Archive metadata rarely says which 10x kit a library used, and when it does it is often wrong.
The reads themselves are reliable: the barcode read's length and its matches against the 10x
whitelists identify the chemistry, and the same evidence shows when a "10x" submission is
something else — scATAC, a feature-barcode library, V(D)J, or not 10x at all. The
[10x chemistry reference](../10x_chemistry_reference.md) is the table this inference implements.

## Per run, then per sample

1. **Per run.** `RENAME10XRUN` samples reads from each run's FASTQs, matches the barcode read
   against every whitelist in `--wl_basedir`, checks the read geometry against the candidate
   layouts, renames the files to Cell Ranger's `<id>_S1_L001_R1_001.fastq.gz` form, and writes a
   `chemistry.json`. A run whose evidence fits no gene-expression layout fails here, with the
   reason in the task's log, and its sample goes no further.
2. **Per sample.** A sample's run reports are collapsed to one chemistry id. When its runs
   disagree, or none has a usable report (every BAM-derived run, for instance), the sample gets
   no forced chemistry and each aligner falls back on its own detection; the log says which.

Every run's call, its confidence, the read lengths and any warnings are in
`results/run_chemistry.tsv`, and per sample in `results/samples.tsv`.

## What each aligner is told

The one per-sample answer is translated separately for each aligner, and deliberately not into
the same set:

| Aligner | Gets | Why |
|---|---|---|
| STARsolo | `starsolo 10x --wl <id>` for any gene-expression id | All `--wl` decides is the whitelist and the barcode/UMI geometry, so a forced id gives what the wrapper's own detection would — without the failure mode where it matches 200k reads, finds no whitelist above its threshold, and aborts on data that matched perfectly well per run |
| Cell Ranger | `--chemistry ARC-v1` for Multiome gene expression only; `auto` otherwise | Cell Ranger detects ARC-v1 itself and then refuses to proceed unless it was named. The `layout_only` ids (`gex_3pv2_or_5pv1v2`, `gex_3pv3_family`), where one whitelist covers several kits, are left to `auto`, which resolves them correctly |

ATAC ids reach neither: `starsolo 10x` is gene expression only.

## When a run is rejected

A rejection is usually right: most turn out to be correctly labelled data that is simply not
single-cell gene expression. [Archive pathologies](../archive-pathologies.md) lists the
submissions that really were broken, and the separate list of correct data the pipeline rejected
late. To check a claim against the FASTQs yourself, the `reprocess-debug` skill's verification
recipes run the inference script directly on a run's work directory.
