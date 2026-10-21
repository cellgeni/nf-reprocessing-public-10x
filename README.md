# nf-reprocessing-public-10x

Fetches public 10x single-cell datasets from GEO, SRA, ENA or ArrayExpress and reprocesses them
uniformly: archive metadata, download (FASTQ, 10x BAM or `.sra`), Cell Ranger FASTQ naming,
chemistry inferred from the reads, and alignment with STARsolo or Cell Ranger, with one run
report across every sample. Runs on the Sanger farm with `-profile sanger`, or anywhere with
Docker or Singularity.

**Documentation: <https://cellgeni.github.io/nf-reprocessing-public-10x/>**

## Quick start

On the Sanger farm:

```bash
module load cellgen/nextflow/26.04.6

nextflow run main.nf -profile sanger \
  --datasets datasets.tsv --default_specie human --starsolo
```

Anywhere else, choose a container engine and pass the whitelists and references, which only
`-profile sanger` provides:

```bash
nextflow run main.nf -profile docker \
  --datasets datasets.tsv --default_specie human --starsolo \
  --wl_basedir whitelists --star_human_reference /path/to/star/GRCh38
```

`datasets.tsv` is tab-separated, one row per series; an empty `sample_id` means the whole series,
so run `--metaonly` first to see what it holds:

```
dataset_id	sample_id
GSE230685	GSM7232572,GSM7232573
PRJEB37166	ERS4605100,ERS4605101
```

Results go to `results/`; start with the run report in `results/multiqc/`. Run
`nextflow run main.nf --help` for the parameters.

## Documentation

| Page | |
|---|---|
| [Getting started](https://cellgeni.github.io/nf-reprocessing-public-10x/getting-started.html) | Requirements, the datasets table, a first run on the farm or elsewhere |
| [Guides](https://cellgeni.github.io/nf-reprocessing-public-10x/guides/) | Profiles and running outside Sanger, resources and queues, species, chemistry inference, reruns and cleanup |
| [Reference](https://cellgeni.github.io/nf-reprocessing-public-10x/reference/) | Parameters, datasets table, outputs, run report, 10x chemistry reference |
| [Troubleshooting](https://cellgeni.github.io/nf-reprocessing-public-10x/troubleshooting.html) | Error messages and their fixes |
| [Knowledge base](https://cellgeni.github.io/nf-reprocessing-public-10x/knowledge-base.html) | Archive pathologies, failure modes, reporting upstream, post-mortems |
| [Developers](https://cellgeni.github.io/nf-reprocessing-public-10x/developers.html) | Layout, architecture, modules, configuration, tests |

The pages are also readable on GitHub in [`docs/`](docs/).

## Repository

| Path | Contents |
|---|---|
| `main.nf`, `workflow/`, `subworkflows/` | The pipeline |
| `nextflow.config`, `configs/` | Parameters, resources by label, per-process settings, the `sanger` profile |
| `nextflow_schema.json` | Parameter schema: validation and `--help` |
| `modules/cellgeni/`, `modules/local/` | Registry and repo-local modules |
| `lib/`, `notebooks/` | STARsolo resource prediction and the analysis behind it |
| `plugins/nf-lineage-keep/` | Work-dir cleanup plugin (`-profile lineage_keep`) |
| `examples/` | An example datasets table and run script |
| `tests/`, `nf-test.config` | nf-test configuration and the LSF test runner |
| `docs/` | Documentation site |

## Changes

See the [changelog](docs/changelog.md).

## Licence

[GPL-3.0](LICENSE)
