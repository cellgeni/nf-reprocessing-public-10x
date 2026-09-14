# cellgeni/wget10x

## Summary

Downloads a single public sequencing file (BAM, SRA, or FASTQ) for a 10x run via `wget`
and renames it to a consistent per-run naming scheme.

The module:

1. Downloads the file at the given URL using `wget`.
2. Renames the downloaded file based on `meta.type`: BAM files are renamed to `${meta.id}.bam`,
   SRA files are renamed to `${meta.id}`. FASTQ files are left as downloaded.
3. Emits the result on the output channel matching its type (`bam`, `sra`, or `fastq`).

## Inputs

| Name | Type | Description |
|---|---|---|
| `meta.id` | string | Run accession used as the renamed output filename (e.g. `SRR12345678`). |
| `meta.type` | string | Expected file type for this run: `BAM`, `SRA`, or `FASTQ`. Determines how the downloaded file is renamed. |
| `link` | string | URL of the file to download for this run. |

## Outputs

| Name | File(s) | Description |
|---|---|---|
| `sra` | `${meta.id}` | Downloaded SRA file, renamed to the run accession. Emitted only when `meta.type` is `SRA`. |
| `fastq` | `*.f*q*` | Downloaded FASTQ file(s). Emitted only when `meta.type` is `FASTQ`. |
| `bam` | `${meta.id}.bam` | Downloaded BAM file, renamed to `${meta.id}.bam`. Emitted only when `meta.type` is `BAM`. |
| `versions` | `versions.yml` | wget version record. |

## Usage

```nextflow
include { WGET10X } from 'cellgeni/wget10x'

WGET10X(
    channel.of([[id: 'SRR12345678', type: 'BAM'], 'https://example.org/SRR12345678.bam'])
)
```

## License

MIT
