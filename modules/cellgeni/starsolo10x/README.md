# cellgeni/starsolo

## Summary

Aligns 10x single-cell RNA-seq FASTQs to a reference genome using STARsolo,
producing per-sample output directories with count matrices and summary statistics.

The module:

1. Stages FASTQs into a per-sample subdirectory under `fastqs/`.
2. Renames the reference index directory to include the species name.
3. Runs the `starsolo 10x` wrapper to perform alignment and cell barcode/UMI counting (no BAM output).

## Inputs

| Name | Type | Description |
|---|---|---|
| `meta.id` | string | Sample accession used as the output directory name (e.g. `SRR12345678`). |
| `meta.wl` | string | Optional. 10x chemistry id (e.g. `gex_3pv3_family`, `gex_multiome_arc_v1`). |
| `fastqs` | files | Gzipped FASTQ files for the sample. |
| `ref_meta.id` | string | Species/reference identifier used to name the genome directory (e.g. `homo_sapiens`). |
| `reference` | directory | STAR genome index directory. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `mapping` | `${meta.id}/` | STARsolo output directory containing count matrices, barcodes, features, and alignment summaries. |
| `versions` | `versions.yml` | STAR and STARsolo wrapper version record. |

## Chemistry (`meta.wl`)

Without `meta.wl` the wrapper detects the chemistry itself: it subsamples reads, matches
barcodes against every known whitelist, and aborts if none of them clears its threshold.

Set `meta.wl` to a chemistry id and the wrapper takes the whitelist and CB/UMI geometry
from that id instead (`starsolo 10x --wl`), skipping barcode matching. The ids are the
ones `infer_10x_run.py` emits, so a value read from its per-run JSON report can be passed
straight through: `gex_3pv1`, `gex_3pv2_or_5pv1v2`, `gex_3pv3_family`, `gex_3pv4_gemx`,
`gex_5pv3_gemx`, `gex_multiome_arc_v1`. ATAC ids are rejected — this module is gene
expression only.

`meta.wl` also implies `--skip-length-checks`, since the same inference that produced the
id has already validated the read lengths; STARsolo then warns instead of aborting on an
unexpected R1/R2 length. UMI-to-R1 reconciliation still applies, and an R1 too short to
hold the barcode is still fatal. Strand specificity is determined by test alignment either
way — it cannot be read off a whitelist.

Requires wrapper v4.3 or later.

## Usage

```nextflow
include { STARSOLO10X } from 'cellgeni/starsolo10x'

STARSOLO10X(
    channel.of([[id: 'SRR12345678', wl: 'gex_3pv3_family'], file('fastqs/')]),
    channel.of([[id: 'homo_sapiens'], file('genome/')])
)
```

## License

MIT
