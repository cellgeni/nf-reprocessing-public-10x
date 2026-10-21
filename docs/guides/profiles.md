---
title: Profiles and running elsewhere
parent: Guides
nav_order: 1
---

# Profiles and running elsewhere
{: .no_toc }

1. TOC
{:toc}

## What a profile is for

The pipeline itself names no executor, queue or container engine: on its own it would run every
task on the local machine, outside any container, and `main.nf` refuses to start that way. A
profile supplies the rest.

| Profile | Sets |
|---|---|
| `sanger` | The Sanger farm: LSF, Singularity with the cellgeni image cache and `/lustre,/nfs` binds, the farm's whitelists and references, and the queues. Self-contained. |
| `singularity` | Singularity, with `autoMounts`. The executor is unchanged. |
| `docker` | Docker, running as your uid. The executor is unchanged. |
| `lineage_keep` | Work-dir cleanup during the run ([Reruns and cleanup](reruns.md)); combine as `-profile sanger,lineage_keep`. |

Profiles combine with commas, and a `-c` config adds anything they do not cover. On the farm,
`-profile sanger` alone is a complete setup; anywhere else, pick an engine and pass the
references.

## How resources become queues

This follows nf-core. A module says how big it is with a label, the pipeline turns the label
into cpus, memory and time, and a site profile picks the queue from the resolved request —
never from the process name. [Resources and queues](resources.md) lists the labels and what each
asks for.

Under `-profile sanger` the queue is chosen like this:

| Request | Queue |
|---|---|
| a task labelled `process_network` (metadata fetch, downloads) | `transfer` |
| memory above 683 GB | `hugemem` |
| time up to 12 h | `normal` |
| time up to 48 h | `long` |
| time up to 7 days, or more | `week`, `basement` |

683 GB is the hard memory limit of `normal` and `long`, which refuse a bigger request at
submission. Only `normal` backfills, which is why STARsolo's run limit is predicted per sample
rather than fixed.

## Running outside the Sanger farm

Pick a container engine and pass the references, which only `sanger` provides.

### Whitelists

Chemistry inference matches each run's barcodes against the 10x whitelists, so it needs a
directory of them under their Cell Ranger names (`737K-august-2016.txt`, `3M-february-2018.txt`,
…). The STARsolo image ships one; copy it out:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/out quay.io/cellgeni/starsolo:v4.4 \
  cp -r /opt/STARsolo/data/whitelists /out/
# or
singularity exec docker://quay.io/cellgeni/starsolo:v4.4 cp -r /opt/STARsolo/data/whitelists .
```

and pass `--wl_basedir whitelists`. STARsolo itself reads its whitelists from inside the image;
only chemistry inference needs the copy.

{: .note }
The image's set lacks `737K-cratac-v1.txt`, the scATAC whitelist. Inference logs a warning and
carries on; an scATAC run is then still rejected, but as matching no supported layout rather than
by name.

### References

| Aligner | Parameters | Needs |
|---|---|---|
| STARsolo (`--starsolo`) | `--star_human_reference`, `--star_mouse_reference` | A STAR `genomeGenerate` directory compatible with STAR 2.7.10a, the version in `quay.io/cellgeni/starsolo:v4.4`. The farm's are built from Cell Ranger's 2020-A GRCh38 and mm10. |
| Cell Ranger (`--cellranger`) | `--cr_human_reference`, `--cr_mouse_reference` | A Cell Ranger reference directory. |

At least one per aligner is required. With only one, samples of the other species are not
aligned, and a warning at launch says so.

### Running it

```bash
nextflow run main.nf -profile docker \
  --datasets datasets.tsv --default_specie human --starsolo \
  --wl_basedir whitelists \
  --star_human_reference /path/to/star/GRCh38 \
  -c local.config
```

With Singularity, set `NXF_SINGULARITY_CACHEDIR` so the images — Cell Ranger's is several GB —
are pulled once rather than into every work directory:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity-cache
nextflow run main.nf -profile singularity ...
```

### Resource caps

The heavy steps ask for 16 cpus and STARsolo for 48 GB or more, and the local executor refuses a
task bigger than the machine. Cap the requests in `local.config`:

```groovy
process.resourceLimits = [cpus: 8, memory: 60.GB, time: 48.h]
```

A human STAR index needs about 32 GB of RAM to load, so below that STARsolo cannot run at all.
When the pipeline runs inside a batch job, also tell the local executor how much of the node it
has:

```groovy
executor.cpus   = 16
executor.memory = 62.GB
```

### A cluster

Give the executor and a queue rule in the same file, choosing the queue from `task.memory` and
`task.time` the way [`configs/sanger.config`](https://github.com/cellgeni/nf-reprocessing-public-10x/blob/main/configs/sanger.config)
does. For SLURM, say:

```groovy
process {
    executor       = 'slurm'
    queue          = { task.time > 24.h ? 'long' : 'short' }
    resourceLimits = [cpus: 64, memory: 500.GB, time: 7.d]
    // Optional: send the metadata and download tasks to a partition with outbound network
    withLabel: 'process_network' { queue = 'transfer' }
}
```

and run with `-profile singularity -c cluster.config`.

{: .note }
Only the tasks labelled `process_network` (FETCH10XMETA and WGET10X) fetch from the archives;
everything after them reads local files. If your compute nodes have no internet access, route
that label to a partition or host that does.
