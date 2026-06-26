# cellgeni/rename10xrun

## Summary

Infers the 10x chemistry and FASTQ layout for a single sequencing run from 2-4 input
FASTQ files, then emits Cell Ranger-style renamed FASTQs and a per-run chemistry JSON
manifest.

The module:

1. Runs `infer_10x_run.py` on the staged FASTQs, which samples up to 200,000 records
   per file and matches barcodes against 10x whitelist files to detect chemistry
   (GEX 3'/5', ATAC, Multiome, GEM-X, etc.).
2. Assigns each FASTQ to a read role (R1, R2, I1, I2) based on barcode geometry and
   read length.
3. Emits renamed `SAMPLE_S1_L001_R1_001.fastq.gz`-style output files and a per-run
   `{run_id}.chemistry.json` manifest recording the detected chemistry, CB/UMI lengths,
   strand hint, and output file mapping.

## Inputs

| Name | Type | Description |
|---|---|---|
| `meta.id` | string | Run accession used as the run ID and output file prefix (e.g. `SRR12345678`). |
| `fastqs` | files | 2-4 gzipped FASTQ files for the run. |
| `whitelist_dir` | directory | Directory containing 10x barcode whitelist files. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `reads` | `*_R*_001.fastq.gz` | Cell Ranger-style renamed biological read FASTQs (R1, R2). |
| `index` | `*_I*_001.fastq.gz` | Cell Ranger-style renamed index read FASTQs (I1, I2), if present. |
| `chemistry` | `${meta.id}.chemistry.json` | Per-run JSON manifest with detected chemistry, modality, CB/UMI lengths, strand hint, and output file mapping. |
| `versions` | `versions.yml` | Python and script version record. |

## Usage

```nextflow
include { RENAME10XRUN } from 'cellgeni/rename10xrun'

RENAME10XRUN(
    channel.of([[id: 'SRR12345678'], [file('SRR12345678_1.fastq.gz'), file('SRR12345678_2.fastq.gz')]]),
    file('whitelists/')
)
```

## License

MIT
