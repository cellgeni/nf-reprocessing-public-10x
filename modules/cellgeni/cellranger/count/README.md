# cellgeni/cellranger/count

## Summary

Runs Cell Ranger count on 10x single-cell RNA-seq FASTQs to align reads,
generate feature-barcode matrices, and produce per-sample output directories.

The module:

1. Stages FASTQs into a per-sample subdirectory under `fastqs/`.
2. Runs `cellranger count` with the provided reference genome directory.
3. Produces a per-sample output directory containing `outs/` with count matrices and summary statistics.

## Inputs

| Name | Type | Description |
|---|---|---|
| `meta.id` | string | Sample accession used as the output directory name and as the FASTQ filename prefix. Must match the prefix of the FASTQ files in `fastqdir` (e.g. `SRR12345678_S1_L001_R1_001.fastq.gz`) so that Cell Ranger can locate the correct files via `--sample`. |
| `fastqdir` | directory | Directory containing FASTQ files for the sample. FASTQ filenames must be prefixed with `meta.id`. |
| `ref_meta.id` | string | Reference identifier (e.g. `human`). |
| `reference` | directory | Cell Ranger reference genome directory. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `mapping` | `${meta.id}/` | Cell Ranger output directory containing `outs/` with count matrices, barcodes, features, and summary statistics. |
| `versions` | `versions.yml` | Cell Ranger version record. |

## Usage

```nextflow
include { CELLRANGER_COUNT } from 'cellgeni/cellranger/count'

CELLRANGER_COUNT(
    channel.of([[id: 'SRR12345678'], file('fastqs/')]),
    channel.of([[id: 'human'], file('reference/')])
)
```

## License

MIT
