# cellgeni/sra2fastq

## Summary

Converts a local SRA file to split gzipped FASTQ files using `parallel-fastq-dump`,
producing one FASTQ per read slot ready for downstream processing.

The module:

1. Runs `parallel-fastq-dump` with `--split-files` to produce per-slot FASTQ files
   (e.g. `SRR12345678_1.fastq`, `SRR12345678_2.fastq`).
2. Compresses each output FASTQ in parallel using `pigz`.

## Inputs

| Name | Type | Description |
|---|---|---|
| `meta.id` | string | SRA run accession (e.g. `SRR12345678`). |
| `sra` | file | Local SRA file for the run. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `fastq` | `*.fastq.gz` | Split gzipped FASTQ files (e.g. `SRR12345678_1.fastq.gz`, `SRR12345678_2.fastq.gz`). |
| `versions` | `versions.yml` | parallel-fastq-dump and pigz version record. |

## Usage

```nextflow
include { SRA2FASTQ } from 'cellgeni/sra2fastq'

SRA2FASTQ(
    channel.of([[id: 'SRR12345678'], file('SRR12345678.sra')])
)
```

## License

MIT
