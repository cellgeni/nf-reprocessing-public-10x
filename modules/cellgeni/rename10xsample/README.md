# cellgeni/rename10xsample

## Summary

Collects per-run Cell Ranger-style FASTQs for a sample and renames them to a consistent
sample-level naming scheme, ready for downstream alignment.

The module runs `rename_fastqs.py` over all staged FASTQs for the sample, producing a
unified set of `SAMPLE_S1_L001_R1_001.fastq.gz`-style output files with consistent
lane numbering across runs.

## Inputs

| Name | Type | Description |
|---|---|---|
| `meta.id` | string | Sample accession used as the output file prefix (e.g. `GSM1234567`). |
| `fastqs` | files | Cell Ranger-style FASTQ files from one or more runs belonging to the sample. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `reads` | `*_R*_001.fastq.gz` | Sample-level renamed biological read FASTQs (R1, R2). |
| `index` | `*_I*_001.fastq.gz` | Sample-level renamed index read FASTQs (I1, I2), if present. |
| `versions` | `versions.yml` | Python version record. |

## Usage

```nextflow
include { RENAME10XSAMPLE } from 'cellgeni/rename10xsample'

RENAME10XSAMPLE(
    channel.of([[id: 'GSM1234567'], [file('SRR1_R1_001.fastq.gz'), file('SRR1_R2_001.fastq.gz'),
                                     file('SRR2_R1_001.fastq.gz'), file('SRR2_R2_001.fastq.gz')]])
)
```

## License

MIT
