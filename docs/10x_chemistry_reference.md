# 10x Genomics chemistry and FASTQ geometry reference

_Last updated: 2026-06-26_

This document is a future-maintenance reference for detecting and validating 10x Genomics FASTQ layouts in public-data reprocessing pipelines. It is written for run-level validation after `sra2fastq`, `bam2fastq`, or direct ENA/FASTQ download, before runs are collected into final samples.

The most important rule is: **do not infer more than the FASTQs can prove**. Some chemistries share barcode whitelists and read geometry. In those cases, the run can often still be safely renamed to Cell Ranger style, but the exact assay/version must be recorded as `layout_only` or resolved from metadata.

---

## 1. Core Cell Ranger FASTQ naming convention

For Cell Ranger-compatible FASTQs, use:

```text
[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
```

or, without a lane token:

```text
[Sample Name]_S1_[Read Type]_001.fastq.gz
```

For `cellranger count`, `cellranger vdj`, and `cellranger multi`, read types are:

```text
I1  optional sample index read
I2  optional sample index read
R1  Read 1
R2  Read 2
```

For `cellranger-atac`, 10x documents two accepted ATAC naming conventions:

```text
I1, R1, R2, R3
```

where `R2` is the i5/barcode read and `R3` is genomic read 2, or:

```text
I1, R1, I2, R2
```

where `I2` is the i5/barcode read and `R2` is genomic read 2.

Pipeline policy:

1. A **run-level script** may temporarily use the run accession as the Cell Ranger sample-name part, for example `SRR123_S1_L001_R1_001.fastq.gz`.
2. A **sample-level script** can later remap run prefixes to final sample prefixes, for example `SAMPLE_A_S1_L002_R1_001.fastq.gz`.
3. STARsolo GEX should only receive the GEX `R1/R2` FASTQs. ATAC, V(D)J-only, Feature Barcode-only, and Flex/probe layouts must not be silently passed as simple GEX.

References: 10x Cell Ranger “Specifying input FASTQ files” and 10x Cell Ranger ATAC “Specifying input FASTQ files”.

---

## 2. Master chemistry table for run-level inference

### 2.1 Simple GEX layouts that can be renamed safely

| Internal chemistry id | 10x family | Modality | Whitelist / barcode evidence | CB length | UMI length | Barcode/UMI location | Typical read lengths seen | STARsolo simple GEX? | Exact-version confidence | Notes |
|---|---|---:|---|---:|---:|---|---|---|---|---|
| `gex_3pv1` | GemCode/Chromium Single Cell 3′ v1 | GEX | `737K-april-2014_rc.txt` | 14 | 5 or 10 measured | Split technical reads: separate CB read + separate UMI read; synthetic STARsolo `R1 = CB+UMI` | cDNA commonly ~98; CB 14; UMI 5 or 10; sample index often 8 | Yes, after CB+UMI reconstruction | Unique at layout level | Early v1 data can have 5 nt UMI. Do not assume 10. If the separate UMI read is missing, fail closed. |
| `gex_3pv2_or_5pv1v2` | Chromium Single Cell 3′ v2 **or** 5′ v1/v2 | GEX | `737K-august-2016.txt` | 16 | 10 | `R1[1:16] = CB`, `R1[17:26] = UMI` | R1 usually 26; biological R2 variable, typically >= 40 | Yes, but strand needs metadata/test alignment | Layout-only | FASTQ geometry and whitelist cannot distinguish 3′ v2 from 5′ v1/v2. |
| `gex_3pv3_family` | Chromium/Next GEM Single Cell 3′ v3/v3.1/LT/HT | GEX | `3M-february-2018.txt` | 16 | 12 | `R1[1:16] = CB`, `R1[17:28] = UMI` | R1 usually 28; biological R2 variable | Yes | Family-level only | v3, v3.1, LT, and HT share the same simple FASTQ geometry. |
| `gex_3pv4_gemx` | GEM-X Universal 3′ v4 | GEX | `3M-3pgex-may-2023.txt` | 16 | 12 | `R1[1:16] = CB`, `R1[17:28] = UMI` | R1 usually 28; biological R2 variable | Yes | Unique enough for current whitelist family | GEM-X v4 has a distinct whitelist; no HT equivalent to 3′ v3 HT. |
| `gex_5pv3_gemx` | GEM-X Universal 5′ v3 | GEX | `3M-5pgex-jan-2023.txt` | 16 | 12 | `R1[1:16] = CB`, `R1[17:28] = UMI` | R1 usually 28; biological R2 variable | Yes, strand reverse | Unique enough for current whitelist family | GEM-X 5′ v3 uses 16 bp CB + 12 bp UMI. |
| `gex_multiome_arc_v1` | Chromium Single Cell Multiome Gene Expression ARC v1 | GEX | `gex_737K-arc-v1.txt` or legacy/shared `737K-arc-v1.txt` | 16 | 12 | `R1[1:16] = CB`, `R1[17:28] = UMI` | Official guide: R1 28, i7 10, i5 10, R2 90 | Yes | Unique only if separated from ATAC | Must not mix with Multiome ATAC. ARC GEX and ATAC share similar barcode-space context but barcode position/read geometry differs. |

### 2.2 ATAC layouts to detect and reject from STARsolo GEX

| Internal chemistry id | 10x family | Modality | Whitelist / barcode evidence | CB length | UMI length | Barcode location | Typical read lengths | STARsolo simple GEX? | Notes |
|---|---|---:|---|---:|---:|---|---|---|---|
| `atac_v1_v2` | Chromium Single Cell ATAC v1/v1.1/v2 | ATAC | `737K-cratac-v1.txt` | 16 | none | i5/Index2 barcode read | R1 50, I1/i7 8, I2/i5 16, R2 50 | No | ATAC genomic reads can be equal length; rely on filenames or explicit roles for R1/R2 orientation if both are 50. |
| `atac_multiome_arc_v1` | Chromium Single Cell Multiome ATAC ARC v1 | ATAC | `atac_737K-arc-v1.txt` or legacy/shared `737K-arc-v1.txt` | 16 | none | i5/Index2 barcode read; often first 16 nt, or nt 9-24 if a 24-cycle i5 includes 8 dark/spacer cycles | Official guide: Read 1N 50, i7 8, i5 24, Read 2N 49. Custom recipe can appear as barcode read length 16 after dark cycles. Public data often appears as R1=50, R2=16/24, R3=49. | No | This is the main “do not mix with GEX” hazard. A 16/24 nt ATAC barcode-index read must never be accepted as GEX R1. |

### 2.3 Supported by Cell Ranger but not safe for simple FASTQ-only GEX inference

| Family | Cell Ranger support status | FASTQ-only inference status | What to do |
|---|---|---|---|
| 5′ V(D)J-only | Supported by `cellranger vdj` / `cellranger multi` depending on version | Not simple STARsolo GEX; may share 5′ barcode/UMI geometry | Require metadata/libraries.csv or explicit library type. Do not feed as GEX count input. |
| 5′ GEX + V(D)J | Supported combinations exist | GEX FASTQs can be processed as 5′ GEX; V(D)J FASTQs need V(D)J pipeline/library type | Separate libraries using metadata; validate each run by library type. |
| 3′ / 5′ Antibody Capture | Supported in combinations and sometimes feature-only modes | Feature Barcode FASTQs have barcode+UMI but distinct feature read/capture structure | Require `libraries.csv` / feature reference; do not infer as standalone GEX without metadata. |
| 3′ / 5′ CRISPR Guide Capture | Supported in combinations | Feature Barcode-like, not ordinary transcriptome GEX | Require feature reference and library type metadata. |
| Antigen Capture | Supported in some 5′ immune profiling combinations | Not simple GEX | Require metadata/libraries.csv. |
| Flex v1 / Flex v2 Gene Expression | Supported by Cell Ranger v10.0 for listed combinations | Probe-based / not a standard simple CB+UMI-in-R1 transcriptome layout for STARsolo | Fail closed in STARsolo-oriented pipelines unless a dedicated Flex handler exists. |
| CellPlex / cell multiplexing / OCM / hashing | Multiplexing solutions supported under specific version/chemistry constraints | Multiplexing tags can share or add feature-barcode libraries | Require library metadata and Cell Ranger multi configuration. |
| Targeted Gene Expression | Support varies by chemistry/version | The FASTQ layout may look like its parent GEX assay, but target-panel metadata is required | Treat as GEX layout for renaming only if metadata confirms targeted GEX; downstream requires appropriate reference/panel setup. |

---

## 3. Detailed chemistry notes

### 3.1 3′ v1: split CB and UMI

3′ v1 is the only GEX layout in this document where CB and UMI can be separate FASTQ reads. The historical/typical structure is:

```text
cDNA biological read: ~98 cycles
cell barcode read:   14 cycles, matches 737K-april-2014_rc.txt
sample index read:    8 cycles
UMI read:             10 cycles in common documentation, but 5 nt occurs in real public v1 data
```

Run-level inference policy:

1. Find the 14 nt CB read by whitelist match against `737K-april-2014_rc.txt`.
2. Find a separate UMI read with length in the allowed set, default `{5, 10}`.
3. Find the long cDNA read.
4. Validate equal read counts and read-ID concordance.
5. Emit synthetic `R1 = CB + UMI`; emit cDNA as `R2`.
6. If the UMI read is not present, fail closed.

STARsolo parameters after reconstruction:

```text
--soloType CB_UMI_Simple
--soloCBstart 1
--soloCBlen 14
--soloUMIstart 15
--soloUMIlen 5 or 10, measured from the data
```

### 3.2 3′ v2 vs 5′ v1/v2 ambiguity

Both use the `737K-august-2016.txt` whitelist, 16 bp CB, and 10 bp UMI. FASTQ geometry alone cannot reliably distinguish:

```text
3′ v2 GEX
5′ v1/v2 GEX
```

Renaming is still safe if one CB+UMI read and one biological read are identified, but exact chemistry/strand should be `layout_only` or decided from metadata or a strand test. For STARsolo, 3′ generally uses `--soloStrand Forward`; 5′ generally uses `--soloStrand Reverse`.

### 3.3 3′ v3/v3.1/LT/HT

These share the same basic FASTQ geometry and `3M-february-2018.txt` whitelist:

```text
CB = 16 bp
UMI = 12 bp
R1 = 28 cycles usually, but can be over-sequenced
R2 = biological read, usually >= 40 bp
```

Use family-level labels, not exact sub-version labels, unless metadata identifies LT/HT/v3/v3.1.

### 3.4 GEM-X Universal 3′ v4

The v4 layout is the same simple CB+UMI structure as 3′ v3/v3.1:

```text
CB = 16 bp
UMI = 12 bp
R1 = 28 cycles usually
```

The key differentiator is the v4 whitelist:

```text
3M-3pgex-may-2023.txt
```

### 3.5 GEM-X Universal 5′ v3

The 5′ v3 layout uses:

```text
CB = 16 bp
UMI = 12 bp
R1 = 28 cycles usually
biological read = R2
```

Use `3M-5pgex-jan-2023.txt` where available. For STARsolo, the strand is expected to be reverse for 5′ gene expression.

### 3.6 Multiome GEX vs Multiome ATAC

Multiome is two libraries from the same cells:

```text
Gene Expression library: simple GEX-like CB+UMI in R1
ATAC library: barcode in index/i5-like read, no UMI, genomic reads in R1/R2 or R1/R3
```

Official Rev F sequencing summary:

```text
Multiome GEX:
  R1  = 28 cycles  (16 bp 10x barcode + 12 bp UMI)
  i7  = 10 cycles
  i5  = 10 cycles
  R2  = 90 cycles  (insert)

Multiome ATAC:
  Read 1N = 50 cycles  (genomic/open chromatin)
  i7      = 8 cycles   (sample index)
  i5      = 24 cycles  (10x barcode + spacer/dark-cycle behavior)
  Read 2N = 49 cycles  (genomic/open chromatin)
```

Public data often appears in either of these ATAC naming forms:

```text
R1=50, R2=16 or 24, R3=49
```

or:

```text
I1=8, R1=50, I2=16 or 24, R2=49
```

Run-level policy:

1. If both GEX and ATAC barcode evidence exists in one run group, fail by default.
2. Only allow `--prefer-gex` if the GEX biological read can be separated unambiguously, for example R2 90/91 while ATAC genomic reads are 49/50.
3. Do not use `--prefer-gex` for Wang-style over-sequenced data where both GEX and ATAC contain 150 bp biological/genomic reads and geometry is insufficient.
4. A 16/24 bp ARC barcode-index read must never be treated as a short GEX R1.

### 3.7 Single Cell ATAC

Standalone 10x scATAC v1/v1.1/v2 uses a different whitelist from GEX:

```text
737K-cratac-v1.txt
```

Typical structure:

```text
R1  = 50 cycles  genomic/open chromatin
I1  = 8 cycles   sample index
I2  = 16 cycles  10x cell barcode
R2  = 50 cycles  genomic/open chromatin
```

Cell Ranger ATAC can also accept BCL-convert-style names where the i5/barcode read is `R2` and genomic read 2 is `R3`. In a STARsolo GEX pipeline, scATAC must be rejected.

### 3.8 Feature Barcode / CRISPR Guide Capture

3′ Feature Barcode introduced additional bead oligos for antibody/CRISPR features. The feature barcode read contains CB+UMI, but it is not the same library type as transcriptome GEX. For 3′ v3/v4 Feature Barcode, the feature library often needs a translation whitelist/mapping file rather than the ordinary GEX whitelist.

Run-level policy:

1. If metadata says `Gene Expression`, process as GEX.
2. If metadata says `Antibody Capture`, `CRISPR Guide Capture`, `Antigen Capture`, or feature-only, do not classify as GEX just because CB+UMI is present.
3. Require `libraries.csv`/feature reference metadata for Cell Ranger `count`/`multi`.

### 3.9 V(D)J

V(D)J libraries belong to `cellranger vdj` or `cellranger multi`, not STARsolo GEX. They can share 5′ barcode/UMI structures with GEX but the biological reads are targeted immune-receptor reads. They must be separated by metadata or by the library construction step; FASTQ read lengths alone are not a safe classifier.

### 3.10 Flex v1/v2

Flex v1/v2 Gene Expression is supported by recent Cell Ranger releases, including combinations with Antibody Capture and/or CRISPR. It is probe-based and should not be handled as a simple STARsolo `CB_UMI_Simple` transcriptome run without a dedicated Flex implementation.

Run-level policy:

```text
If suspected Flex and no dedicated handler exists: fail closed.
```

---

## 4. Whitelist files and intended usage

| Whitelist filename | Intended family | CB length | Notes |
|---|---|---:|---|
| `737K-april-2014_rc.txt` | 3′ v1 | 14 | Used for split CB read. UMI read is separate and must be measured; allow 5 or 10 nt. |
| `737K-august-2016.txt` | 3′ v2 and 5′ v1/v2 | 16 | Ambiguous between 3′ v2 and 5′ v1/v2. |
| `3M-february-2018.txt` | 3′ v3/v3.1/LT/HT | 16 | Family-level only from FASTQs. |
| `3M-3pgex-may-2023.txt` | GEM-X Universal 3′ v4 | 16 | v4 GEX barcode whitelist. |
| `3M-5pgex-jan-2023.txt` | GEM-X Universal 5′ v3 | 16 | 5′ v3 GEX barcode whitelist. |
| `gex_737K-arc-v1.txt` | Multiome GEX ARC v1 | 16 | Prefer this over shared `737K-arc-v1.txt` when present. |
| `atac_737K-arc-v1.txt` | Multiome ATAC ARC v1 | 16 | Barcode in i5/index-like read, not GEX R1. |
| `737K-arc-v1.txt` | Legacy/shared ARC whitelist | 16 | Dangerous unless geometry separates GEX from ATAC. |
| `737K-cratac-v1.txt` | Single Cell ATAC v1/v1.1/v2 | 16 | Different from GEX barcodes; ATAC-only, no UMI. |
| `translation_3M-february-2018.txt.gz` | 3′ Feature Barcode translation | 16 | Feature barcode mapping; do not confuse with ordinary GEX whitelist. |
| `translation_3M-3pgex-may-2023.txt.gz` | 3′ v4 Feature Barcode translation | 16 | Feature barcode mapping for v4 family. |

---

## 5. Read-length signatures to validate across collected runs

### 5.1 GEX-safe signatures

| Signature | Interpretation | Allowed grouping behavior |
|---|---|---|
| `R1=14, UMI=5/10, R2>=40` before reconstruction | 3′ v1 split CB/UMI | Run-level script reconstructs `R1=19` or `R1=24`; validator should compare post-renaming output lengths. |
| `R1=26, R2>=40` | 3′ v2 or 5′ v1/v2 | Group only if chemistry/layout group is identical. Strand unresolved unless metadata/test alignment exists. |
| `R1=28, R2>=40` | 3′ v3/v4, 5′ v3, or Multiome GEX depending on whitelist | Group by whitelist family. |
| `R1>28` with whitelist at offset 0 and full CB+UMI present | Over-sequenced GEX R1 | Safe for barcode extraction if no conflicting modality; record actual R1 length. |
| `I1=10, I2=10, R1=28, R2=90/91` | Multiome GEX | GEX-only. Do not mix with ATAC. |

### 5.2 ATAC signatures to reject from GEX

| Signature | Interpretation | GEX action |
|---|---|---|
| `R1=50, R2=16, R3=49` | Multiome ATAC, barcode in `R2` | Reject from STARsolo GEX; optionally rename for Cell Ranger ARC/ATAC-specific path. |
| `R1=50, R2=24, R3=49` | Multiome ATAC with 24-cycle i5/barcode read | Reject from STARsolo GEX. |
| `I1=8, R1=50, I2=16, R2=49/50` | Multiome or standalone ATAC | Reject from STARsolo GEX. |
| `I1=8, R1=50, I2=24, R2=49` | Multiome ATAC official/custom recipe data | Reject from STARsolo GEX. |
| `I1=24, R1=150, R2=8, R3=150` | Wang-style ATAC deposit | Reject from STARsolo GEX. |

---

## 6. Observed public Multiome read-length examples

These examples came from `cellgeni/multiome-processing/data/readlengths` and are useful regression cases for run-level validators.

| Dataset/readlength source | Library | Observed columns/read lengths | Expected action |
|---|---|---|---|
| `Artyunyan2023` | GEX | `I1=10, I2=10, R1=28, R2=90` | Accept as Multiome GEX. |
| `Artyunyan2023` | ATAC | `I1=8, R1=50, R2=16, R3=49` | Reject from GEX; barcode read is ATAC R2. |
| `Ma2022` | GEX | `I1=10, I2=10, R1=28, R2=90` | Accept as Multiome GEX. |
| `Ma2022` | ATAC | `I1=8, I2=24, R1=50, R2=49` | Reject from GEX; ATAC barcode is I2. |
| `Wang2022` | GEX | `I1=10, I2=10, R1=150, R2=150` | Accept as over-sequenced GEX only if isolated from ATAC and whitelist confirms R1. |
| `Wang2022` | ATAC | `I1=24, R1=150, R2=8, R3=150` | Reject from GEX. Mixed Wang-style GEX+ATAC should fail closed because both include 150 bp long reads. |
| `Kazu2023` | GEX | `I1=10, I2=10, R1=28, R2=90` | Accept as GEX. |
| `Kazu2023` | GEX | `I1=10, I2=10, R1=28, R2=81` | Accept as GEX, but sample-level validation should require same R2 length unless an explicit delta is allowed. |
| `Kazu2023` | ATAC | `I1=8, R1=50, R2=16, R3=50` | Reject from GEX. |
| `Kazu2023` | ATAC | `I1=8, R1=50, R2=24, R3=49` | Reject from GEX. |
| `GSE268630` | ATAC | generic columns `1=8, 2=50, 3=16, 4=49` | Reject from GEX; infer roles by geometry if names are generic. |
| `GSE268630` | GEX | generic columns `1=28, 2=151` | Accept as GEX when whitelist confirms the 28 nt barcode read. |
| `GSE268630` | ATAC | generic columns `1=8, 2=50, 3=24, 4=49` | Reject from GEX. |
| `MassoniBadosa2024` | GEX/RNA | `R1=28, R2=90` | Accept as GEX. |
| `MassoniBadosa2024` | ATAC | `R1=50, R2=16, R3=49` | Reject from GEX. |
| `MassoniBadosa2024` | ATAC | `R1=50, R2=24, R3=49` | Reject from GEX. |

---

## 7. Validation checklist for a run-level script

A run-level script should be invoked on exactly one run/mate set: usually 2-4 FASTQ files. It should fail if a sample-level mixture has already been created.

Minimum hard checks:

1. Input count is 2-4 FASTQs.
2. Compression is detected by magic bytes, not extension.
3. FASTQ records are structurally valid: header starts with `@`, plus line starts with `+`, sequence and quality lengths match, no truncated records.
4. Read lengths are stable in the sampled reads.
5. Barcode candidate is identified by whitelist content, not filename alone.
6. Non-v1 GEX barcode read must contain full `CB+UMI`.
7. v1 must have CB read and UMI read; v1 UMI length must be measured.
8. All emitted read roles have equal exact read counts.
9. First N read IDs match across emitted mates, after normalizing common `/1`, `/2`, `.1`, `.2` suffixes.
10. No long unassigned FASTQ is ignored.
11. Mixed GEX+ATAC evidence fails unless an explicit and safe selection mode is requested.
12. ATAC-only runs fail under `--require-modality gex` or `--require-gex`.

Suggested JSON fields if a run-level report is emitted:

```json
{
  "run_id": "SRR123",
  "ok": true,
  "modality": "gex",
  "chemistry_id": "gex_3pv3_family",
  "chemistry_group": "gex_3pv3_family",
  "layout_confidence": "layout_only|unique",
  "cb_len": 16,
  "umi_len": 12,
  "read_lengths": {"R1": 28, "R2": 91},
  "outputs": ["SRR123_S1_L001_R1_001.fastq.gz", "SRR123_S1_L001_R2_001.fastq.gz"],
  "warnings": []
}
```

A sample-level validator should not require JSON reports from other pipelines. It should be able to validate already-renamed Cell Ranger FASTQs directly by parsing names and sampling read lengths. Optional chemistry checking can be enabled if whitelist files are available.

---

## 8. Known non-identifiability and fail-closed cases

| Case | Why it is ambiguous | Safe behavior |
|---|---|---|
| 3′ v2 vs 5′ v1/v2 | Same `737K-august-2016` whitelist and 16+10 CB/UMI layout | Rename as layout family; require metadata/test alignment for exact chemistry/strand. |
| 3′ v3 vs v3.1 vs LT vs HT | Same whitelist and simple geometry | Rename as family; use metadata for exact label. |
| 5′ GEX vs 5′ V(D)J | Can share 5′ barcode/UMI geometry but different biological target | Require library type metadata. |
| GEX + Feature Barcode | Feature libraries can contain CB+UMI but are not transcriptome reads | Require libraries.csv/feature reference. |
| Flex | Probe-based Cell Ranger chemistry; not simple STARsolo GEX | Dedicated Flex handler or fail closed. |
| Mixed Multiome GEX+ATAC | ARC barcode evidence can be present in both libraries; ATAC barcode is index-like | Fail by default; only select GEX when read lengths clearly separate GEX R2 from ATAC genomic reads. |
| Wang-style over-sequenced Multiome | GEX and ATAC can both have 150 bp long reads | Do not use read length alone; fail closed unless metadata/lane separation is explicit. |

---

## 9. Recommended pipeline policy

### Run-level step

Run on one run only:

```bash
infer_10x_run.py \
  --fastqs <2-4 FASTQs from one run> \
  --run-id <RUN_ID> \
  --whitelist-dir <10x whitelist dir> \
  --outdir <run-fastq-dir> \
  --require-modality gex
```

Temporary output names can use the run accession:

```text
SRR123_S1_L001_R1_001.fastq.gz
SRR123_S1_L001_R2_001.fastq.gz
```

### Sample-level step

Validate and collect already-renamed run FASTQs:

```bash
validate_10x_runs.py \
  --fastqs <all run-level CR-style FASTQs for one sample> \
  --sample-id <FINAL_SAMPLE_ID> \
  --outdir <sample-fastq-dir> \
  --require-modality gex \
  --ignore-index-reads
```

Default compatibility check can be read-length only. Optional whitelist-based chemistry checking can be added for internal pipelines, but should not be required for reuse by other pipelines that already have valid CR-style FASTQs.

---

## 10. Source links

Primary/current documentation:

- 10x Genomics Cell Ranger supported libraries: https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/supported-libraries
- 10x Genomics Cell Ranger FASTQ naming: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-specifying-fastqs
- 10x Genomics Cell Ranger ATAC FASTQ naming: https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/inputs/specifying-input-fastq-files
- 10x Multiome Rev F user guide: https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevF.pdf

Library-structure references:

- 10x 3′ v1: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html
- 10x 3′ v2/v3/v3.1/v4: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html
- 10x 5′ GEX: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html
- 10x 3′ Feature Barcode: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html
- 10x 5′ V(D)J/Feature Barcode: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5vdjfb.html
- 10x Single Cell ATAC: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html
- 10x Multiome ATAC + GEX: https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html

Internal/pipeline references:

- `infer_10x_run.py` chemistry table and inference policy.
- `validate_10x_runs.py` sample-level compatibility checks.
- `cellgeni/multiome-processing/data/readlengths` observed public Multiome read-length regression cases.
- `sra_to_10x_fastq_gz.sh` historical v1 CB+UMI reconstruction logic.
- `platform_10x.sh` historical STARsolo whitelist and read-length detection logic.
