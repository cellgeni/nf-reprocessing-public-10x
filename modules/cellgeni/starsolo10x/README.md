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
| `fastqs` | files | Gzipped FASTQ files for the sample. |
| `ref_meta.id` | string | Species/reference identifier used to name the genome directory (e.g. `homo_sapiens`). |
| `reference` | directory | STAR genome index directory. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `mapping` | `${meta.id}/` | STARsolo output directory containing count matrices, barcodes, features, and alignment summaries. |
| `versions` | `versions.yml` | STAR and STARsolo wrapper version record. |

## Usage

```nextflow
include { STARSOLO10X } from 'cellgeni/starsolo10x'

STARSOLO10X(
    channel.of([[id: 'SRR12345678'], file('fastqs/')]),
    channel.of([[id: 'homo_sapiens'], file('genome/')])
)
```

## License

MIT
