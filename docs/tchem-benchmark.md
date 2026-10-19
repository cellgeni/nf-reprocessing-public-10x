# tchem on our runs

[TusslChemistry](https://github.com/Tussl-Co/TusslChemistry) (`tchem`) is a Rust tool, made
public on 2026-10-14, that names a single-cell library's kit, barcode list, read geometry, cDNA
strand and **library type** from the first 100,000 reads of each file. Its README compares it with
"Sanger's `infer_10x_chemistry.py`", which is a script from this repository.

This page records what tchem 0.1.0 does on runs whose outcome we already know, so the team can
decide whether to take any of it into the pipeline. It compares and does not recommend; the
options are at the end.

Measured 2026-10-19 on 1,227 runs (1,236 run/input combinations). The raw evidence, scripts and
per-run reports are archived at `/nfs/cellgeni/reprocessing-runs/benchmarks/tchem-0.1.0/`
(see its `manifest.txt`).

## Summary

| Question | What we found |
|---|---|
| Does tchem pick the same barcode list as `RENAME10XRUN`? | Yes, on all 1,042 routine runs it scored, with no disagreement. 967 of them are HIGH. |
| Does it tell 3' from 5' on the shared 737K list? | On 157 of 209 such runs (75%). It agrees with STARsolo's test alignment on 156; the one disagreement is a CITE-seq library. The other 52 runs it leaves ambiguous, and STARsolo resolves them. |
| Does it spot V(D)J libraries? | All 36 labelled V(D)J runs, and 41 of 42 routine runs whose ENA library name says TCR or VDJ. |
| Does it spot feature-barcode libraries? | On 3' kits it finds 17 of 23 labelled ones. On **5' kits (CSP, ADT, CITE, hashtag) it finds 0 of 23**: it calls them gene expression. |
| Does it call any real gene-expression run non-GEX? | No: 61 of 61 labelled GEX runs, and 244 of 244 routine runs named as GEX, are GEX. |
| What does it see that our pipeline does not? | **25 of 222 samples** in two recent runs contain V(D)J or tag runs. 22 of these samples are GEX with non-GEX runs merged in under one BioSample, and 3 are V(D)J-only samples published as GEX. Per-sample QC cannot see the merged ones. |
| What does it handle worse than we do? | It refuses runs with a 16 bp index read (8 runs) or an ENA orphan third file (1 run). It makes calls on files of 4-176 reads, 13 of them at HIGH. It has no 3' v1, no ATAC and no bz2. |
| Cost | Median 4.3 s CPU and 145 MB RSS per run, plus a 400 MB barcode-list directory in its own layout. |

## What tchem is, next to `RENAME10XRUN`

| | tchem 0.1.0 | `RENAME10XRUN` (`infer_10x_run_recommended.py`) |
|---|---|---|
| Output | A JSON report; optionally an aligner command line. **Never writes FASTQs.** | Read roles, mate-set checks, Cell Ranger-named FASTQs and `chemistry.json` |
| Kits | 7 10x RNA layouts (3' v2/v3/v4, 5' v1-v2/v3, Multiome GEX, Flex), 18 non-10x layouts, 10x on long reads | 8 10x ids, including 3' v1 and two ATAC ids (ATAC is rejected) |
| Barcode match | Exact or one substitution; one `N` allowed | Exact only |
| Sampling | First 100,000 reads of each file | Seeded reservoir of 200,000 from the first 2,000,000 |
| 3' v2 vs 5' v1/v2 | polyT/TSO tail, then reference-free strand: polyA/T, human and mouse chrM, 40 human nuclear genes | Not separated (`layout_only`); STARsolo's test alignment decides |
| Library type | V(D)J (human constant regions), CRISPR scaffold, tag/amplicon (read-start diversity). Reported, not refused. | GEX vs ATAC only |
| Index reads | By content: at most 12 bp, low diversity | By modal length in `{6,7,8,9,10,14,16,24}` |
| Inputs | 1-4 files, of which 1-2 non-index; gzip, plain or BAM | 2-4 files; gzip, bz2 or plain; drops ENA orphan files |

The chemistry table our inference implements is
[10x_chemistry_reference.md](10x_chemistry_reference.md). tchem's is its README's "Supported
layouts" and "How tchem decides".

## What their published comparison tested

tchem's README has a section "Comparison with Sanger's tool". It ran
`infer_10x_chemistry.py` at commit `a861162` (2026-06-26) of this repository, with
`--dry-run --allow-missing-whitelists` (`bench/head_to_head/scripts/select.py`).

**That script is not what the pipeline runs.** `RENAME10XRUN` calls
`infer_10x_run_recommended.py`. `infer_10x_chemistry.py` is an older, stricter design still
sitting in `resources/usr/bin`:
- it samples the first 200,000 reads instead of a reservoir;
- it needs constant read lengths;
- it has no `--prefer-gex`.

Of its 34 errors in their test, about 12 ("variable read lengths") and 3 ("use --prefer-gex")
are cases the live script handles.

What their numbers do show about our approach:

- **Barcode lists:** on 445 experiments our script named a different list from scBaseCount's
  in **0** cases. It made no call on 19.
- **Second opinion:** in 914 runs where tchem disagreed with the archive's label, our script
  sided with tchem in 851 and with the label in none.
- **3' vs 5':** their headline gap (tchem 87%, ours 51%) is by design. Our script returns
  `gex_3pv2_or_5pv1v2` and leaves the strand to STARsolo.

## How we measured

**Two sets of runs, all declared human.**

- **Set A, labelled pathologies:** 170 runs from 25 datasets whose library type or defect is
  recorded in [archive-pathologies.md](archive-pathologies.md) (batches 9-19).
  - Up to 4 samples per library class per dataset, so a dataset's GEX samples serve as
    controls for its V(D)J or feature samples.
  - **Label:** the GEO title (`…_TCR`, `… CSP`, `… GEX`). Where the title says nothing, the
    published matrix decides: `exon_u` < 0.05 with ≤ 20 median features means a feature
    library; `exon_u` ≥ 0.2 with ≥ 500 features means GEX.
- **Set B, routine runs:** every run of the two recent runs that kept `run_chemistry.tsv`:
  - batch 25's velocyto re-run table;
  - `batch_velocyto_rerun`.

  That is 1,057 runs, 222 samples and 51 datasets, each with our `chemistry_id` and the
  sample's published QC.

**Input.** For each run tchem got the first 100,000 reads of every file. They came from the same
place our pipeline downloaded from: the ENA FASTQ or the submitter's original FASTQ, cut after
100,000 records. For runs we took from SRA, the input was `fastq-dump -X 100000 --split-files`,
technical reads included. For the two ENA-defect datasets both routes were run.

**References.**
- **3' vs 5':** the STARsolo wrapper's own decision in each sample's `strand.txt`, i.e. its
  Forward/Reverse test alignment on GeneFull.
  - Not the `Strand` column of `mapping_qc_stats.tsv`; that column is wrong for paired-end 5'
    data (see [Side findings](#side-findings-about-our-own-pipeline)).
  - Samples on a 3'-only whitelist, where the wrapper forces Forward without a test, are
    counted separately.
- **Library type:** the label above, plus ENA's `library_name` where it names the library
  (`…_TCR`, `…VDJ…`, `…_CiteSeq`, `…_ADT`, `…_CM`).

**Not measured:**
- whether a call within a list family (3' v3 vs v3.1) is right;
- mouse or any other species;
- what the merged non-GEX runs did to the affected matrices.

## Results

### Barcode list

Set B, 1,057 runs: 1,042 got a list call and 15 did not (see
[Where tchem refused or failed](#where-tchem-refused-or-failed)).

| Our `chemistry_id` | Same list | Different list | tchem confidence |
|---|---|---|---|
| `gex_3pv3_family` | 747 | 0 | 728 HIGH, 1 MEDIUM, 18 LOW |
| `gex_3pv2_or_5pv1v2` | 209 | 0 | 153 HIGH, 4 MEDIUM, 52 ambiguous (3' v2 vs 5' v1/v2) |
| `gex_multiome_arc_v1` | 46 | 0 | 46 HIGH |
| `gex_5pv3_gemx` | 31 | 0 | 31 HIGH |
| `gex_3pv4_gemx` | 9 | 0 | 9 HIGH |

All 18 LOW calls are lane files of 4-15 reads in PRJEB52216, where a "non-random UMI position"
is an artefact of the sample size.

**tchem's head-only sampling did not cost a list call on any run.** The <10% head-of-file
scores that made us move to a reservoir came from cycle-1 `N`s under exact matching. tchem
tolerates one mismatch and one `N`, which appears to absorb this.

### 3' vs 5'

Set B, runs where STARsolo ran a test alignment:

| Our `chemistry_id` | STARsolo | tchem agrees | tchem disagrees | tchem ambiguous |
|---|---|---|---|---|
| `gex_3pv2_or_5pv1v2` | 3' (44) | 42 | 0 | 2 |
| `gex_3pv2_or_5pv1v2` | 5' (165) | 114 | 1 | 50 |
| `gex_3pv3_family` | 3' | 635 | 0 | 0 |
| `gex_multiome_arc_v1` | 3' | 46 | 0 | 0 |
| `gex_5pv3_gemx` | 5' | 31 | 0 | 0 |
| `gex_3pv4_gemx` | 3' | 9 | 0 | 0 |

- **Forced Forward:** on the 96 runs where the wrapper forced Forward with no test, tchem also
  says 3' on all 96.
- **Per sample:** 139 agree, 4 have no tchem 3'/5' call on any run, 1 is split between its runs,
  and none disagree.
- **The one disagreement:** `SRR27765505`, ENA library name `…_P9R_CiteSeq`.
  - It is an antibody library inside a 5' sample. tchem has no chrM or nuclear reads to go on,
    so it called 3' v2 from the polyA signal alone, at MEDIUM.
  - STARsolo's test said Reverse, but only 6% vs 30%.
- **Ambiguous runs (52):**
  - 28 have a 26-28 bp barcode read, too short to show a TSO or polyT. tchem's README lists
    this as a known gap.
  - 24 are the V(D)J runs of PRJNA898288.
  - STARsolo resolved 50 of the 52 as 5', with margins no narrower than 4% vs 24%.
  - The other 2 had no signal for either tool.

For 3'/5', tchem adds an alignment-free second opinion that agrees with what we already do. It
does not replace the test alignment, which resolves a quarter more of the 737K runs.

### Library type: labelled datasets (Set A)

| Label | Runs | tchem says |
|---|---|---|
| V(D)J (title `TCR`/`BCR`/`VDJ`) | 36 | **V(D)J 36** |
| Feature barcode, 3' kit | 23 | tag 16, CRISPR 1, **GEX 6** |
| Feature barcode, 5' kit | 23 | **GEX 23** |
| Gene expression | 61 | GEX 61 |
| Unidentified, near-empty matrix (GSE214914, GSE218316) | 8 | **tag 8** |
| Not 10x (GSE222431, bulk Ig 5'RACE) | 4 | no match 4 |

Our pipeline aligned and published every one of these, apart from 2 feature runs with a 20 bp
R2.

**The 5' blind spot.** tchem calls every 5'-kit feature library gene expression:
- GSE216914 CSP;
- GSE219098 ADT;
- GSE221776 CITE;
- GSE229626 CSP;
- GSE231370 Hashtag;
- GSE233046 ADT.

Its tag test looks for a few sequences covering the **start** of the cDNA read. In these
libraries the antibody barcode is not at the start: the top read starts cover only 1.5-6% of
reads.

The six 3' misses:
- **GSE229733 `prot`:** 4 runs. Read starts cover 1.6%, so this looks like the same layout.
- **GSE215253:** the two 20 bp-R2 CRISPR libraries, which tchem calls **Drop-seq at MEDIUM**.
  This is the only wrong layout name in the benchmark. It is not a HIGH call. We rejected these
  two runs.

All 29 missed feature runs share one property: tchem's own `housekeeping_share` is exactly 0,
where the labelled GEX runs score 0.05-0.18. A housekeeping floor alone would still misfire,
though. PRJEB67784, a human GEX dataset in Set B with healthy matrices, has runs scoring
0.0001-0.0005.

**The two "unidentified" datasets** in archive-pathologies.md are tag or amplicon libraries
according to tchem:
- The top 200 read starts cover 96-99% of reads.
- There are no housekeeping or constant-region reads.
- The calls are HIGH on the 3' v3 or v2 list.

This fits their 0.04-0.76% mapping.

### Library type: routine runs (Set B)

tchem types runs as V(D)J or tag in **25 of 222 samples (11%)**, all published as GEX:

| Dataset | Samples | Runs per sample (GEX + other) | ENA library names of the other runs |
|---|---|---|---|
| PRJNA898288 | 3 | 16 GEX + 8 V(D)J | `CLS063_pop5_VDJ_S7_L001`, … |
| PRJNA744851 | 6 | 2-4 GEX + 1 V(D)J | `2169460_SJTRC_InfVax_1_TCR`, … |
| PRJNA944134 | 4 | 1-4 GEX + 1 V(D)J | `PC1tcr` … `PC4tcr` |
| PRJNA908079 | 3 | 1 GEX + 1 V(D)J | `7`, `9`, `11` (library_selection PCR) |
| PRJNA1070700 | 1 | 5 GEX + 3 V(D)J | `…_P9R_TCR`, `…_PepStim_TCR` |
| PRJNA1109397 | 2 | V(D)J only (2) | `Donor1_Blood_TCR`, `Donor1_Tonsil_TCR` |
| GSE181543 | 1 | V(D)J only (1) | GEO title `PBMC2 [VDJ]` |
| PRJEB52216 | 2 | 31-34 GEX + 12-15 tag | `SIGAA12_p`, `SIGAC10_p` |
| PRJNA934842 | 3 | 1 GEX + 1 tag | `S42053_765_CM`, … (cell multiplexing) |

**Mechanism.** Each donor's GEX, V(D)J and feature libraries are separate runs of **one
BioSample**. We group runs by BioSample, so they become one sample and one matrix. This is the
shared-BioSample problem already recorded for GEO datasets (GitHub #6), here in BioProject
datasets.

**Why QC misses it.** The mixed samples' matrices look healthy (`exon_u` 0.18-0.69, 538-5,354
median features). Neither the near-empty screen nor a mapping gate would catch them.

**Checked against ENA's own library names:**
- 41 of the 42 runs named TCR or VDJ are typed V(D)J. The exception, `SRR23906943` `PC5tcr`,
  has a constant-region share of 0.096 against 0.051 housekeeping, under tchem's 3× rule.
- 16 runs named CITE-seq or ADT are typed GEX: the 5' blind spot again (PRJNA744851,
  PRJNA1070700, PRJNA1313799, PRJEB60880).
- No run named as GEX was typed anything else.

### Defects and non-10x data (Set A)

| Dataset | Defect | `RENAME10XRUN` | tchem on ENA's files | tchem on SRA's copy |
|---|---|---|---|---|
| GSE228428 (5 runs) | ENA's export drops the barcode read | rejected: no whitelist hits | no match; hint names the missing barcode read and `fasterq-dump --include-technical` | Multiome GEX, HIGH (5/5) |
| GSE241998 (3 runs) | same | rejected | no match, same hint (3/3) | 3' v3, HIGH (3/3) |
| GSE233338 (4 runs) | index reads only | rejected | refused, exit 1: "0 read files left" | not run |
| GSE222431 (4 runs) | bulk Ig 5'RACE, not 10x | **aligned and published** (97 samples) | no match (4/4) | not run |
| GSE235787 `GSM8086562` | 10x Flex | rejected: no Flex whitelist | 10x Flex v1, HIGH | not run |

GSE235787's 4 sampled Multiome GEX siblings of the Flex sample were all Multiome GEX at HIGH.

### Where tchem refused or failed

On runs our pipeline aligned:

| What | Runs | Why |
|---|---|---|
| 16 bp index read (GSE274751, reads 8/28/16/98) | 8 | Only reads ≤ 12 bp count as index. Three files are left and tchem exits 1. |
| ENA orphan file beside `_1`/`_2` (PRJNA816507 `SRR18365379`, 180 reads) | 1 | Three read files, exit 1. We drop the orphan (issue #10's code). |
| `fastq-dump` "row NULL" (PRJNA524398, BGISEQ) | 6 | Our harness, not tchem. This is the NULL-`NAME` case that `REPROCESS10X_SRAFIX` exists for. |
| Calls on tiny files: 33 runs with < 1,000 reads | 33 | 13 called HIGH (lowest: 4 reads). There is no minimum-read floor. |

The 33 tiny-file runs got calls rather than refusals, so they are not in the 15 runs without a
list call.

From tchem's own scope: no 10x 3' v1, no ATAC, no bz2. All 23 non-Flex barcode lists must be
present, in its own file names. Our `/nfs/cellgeni/STAR/whitelists` cannot be used as it is.

### Cost

| | Median | p95 | Max |
|---|---|---|---|
| tchem CPU per run (user + sys) | 4.3 s | 5.5 s | 6.2 s |
| tchem peak RSS | 145 MB | | 178 MB |
| Fetching the first 100,000 reads of each file | 3 s | 12 s | 291 s |

tchem's README quotes 1.3 s per run, measured on cached files. Ours includes index files and
uncompressed input. `RENAME10XRUN` takes about 40 s CPU per run, most of it a full-file record
count.

Its barcode-list directory is 403 MB including 77 MB of `.packed` files. The lists come from
`scripts/fetch_whitelists.sh`, which checks each one's sha256. The Flex list can only be copied
from a Cell Ranger install.

## Side findings about our own pipeline

- **Merged non-GEX runs** (above): 22 samples in two recent runs carry V(D)J or tag runs inside
  a GEX matrix, and 3 V(D)J-only samples were published as GEX.
- **`mapping_qc_stats.tsv`'s `Strand` and `Paired` columns are wrong for paired-end 5' samples.**
  - All 9 5' GEM-X samples of batch 25's velocyto re-run read `Single` / `Forward` there.
  - The wrapper's `strand.txt` says paired-end mode with Reverse (e.g. 12% vs 65%).
  - `REPROCESS10X_MAPPINGQC` runs `starsolo qc` from the `starsolo:v4.1` image, while the
    alignment uses v4.4.
  - Anything that reads 3'/5' from the run report is affected.
- **GSE214914 and GSE218316**, recorded as `unidentified-library`, read as tag or amplicon
  libraries.

## Options

These are the choices the evidence bears on.

| Option | What it would take | What it buys | What it does not fix |
|---|---|---|---|
| **A. tchem as an advisory step**: one process per run on `RENAME10XRUN`'s output, whose `library.kind`, `cdna_strand` and verdict go into `run_chemistry.tsv` and the run report | A container for the static binary; a second whitelist directory (403 MB, its file layout, Flex from Cell Ranger); about 5 s CPU per run; a new process for `lineage_keep` to classify | V(D)J and 3' tag runs flagged per run, including merged ones; a strand cross-check; a hint for ENA exports missing the barcode read | The 5' feature blind spot; nothing stops the merge unless the pipeline acts on the flag; tchem 0.1.0 is 5 days old, with no conda or container release |
| **B. Port the signals into our Python**: constant-region share (V(D)J), read-start coverage (tag), sgRNA scaffold (CRISPR), possibly the reference-free strand | Its references (GENCODE-derived, Apache-2.0, with attribution); a few hundred lines; our own thresholds checked against Set A and Set B, which are now archived for that | The same flags with no new dependency; room to add what tchem lacks, e.g. barcode-at-offset-10 for 5' feature libraries | Needs its own validation; we would own the maintenance |
| **C. Leave chemistry inference as it is** | Nothing | Nothing new | Merged and V(D)J-only samples keep reaching published matrices; 5' feature libraries stay caught only by the near-empty screen, after alignment |

Independent of the choice: splitting runs by library type before grouping by BioSample is a
`download10x` change that any per-run flag would feed. It would also need its own policy: drop
the non-GEX runs, or publish them separately.

## Reproducing

Everything is under `/nfs/cellgeni/reprocessing-runs/benchmarks/tchem-0.1.0/`:

| File | Contents |
|---|---|
| `sets.tsv` | The two run sets: our call, STARsolo's `strand.txt`, QC, ENA metadata |
| `calls.tsv` | `sets.tsv` with tchem's report fields (`t_*`) joined |
| `summary.md` | Every table `score.py` produced, including each Set A run |
| `reports_json.tar.gz` | tchem's JSON report per run |
| `scripts/` | `build_sets.py` builds the sets; `run_one.sh`/`run_all.sh` fetch and run, through `farm_run.sh` on the normal queue; `score.py` scores |
| `tchem-v0.1.0-….sha256` | The release checksum |
