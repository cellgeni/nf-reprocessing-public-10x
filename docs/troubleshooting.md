---
title: Troubleshooting
nav_order: 5
---

# Troubleshooting
{: .no_toc }

Find the message you got, then the fix. For what a task's exit code means, see
[Failure modes](failure-modes.md).

1. TOC
{:toc}

## The run stops at launch

### `No container engine: run with -profile sanger, singularity or docker`

No profile was given, or none that enables an engine. Add `-profile sanger` on the Sanger farm,
or `-profile singularity` / `-profile docker` elsewhere. Pass it on every `-resume` too; profiles
are not remembered between runs. See [Profiles and running elsewhere](guides/profiles.md).

### `--wl_basedir is required to infer chemistry`

Outside `-profile sanger` there is no default whitelist directory. Copy the one the STARsolo image
ships and pass it — see [Whitelists](guides/profiles.md#whitelists). Not needed with
`--metaonly`.

### `--starsolo needs --star_human_reference and/or --star_mouse_reference`

The same for `--cellranger` with `--cr_human_reference` / `--cr_mouse_reference`. Pass at least
one reference for the aligner you asked for — see [References](guides/profiles.md#references).
With only one, the run goes ahead and a warning says the other species will not be aligned.

### `When --no_infer_specie is set, --default_specie must be provided.`

`--no_infer_specie` gives every sample `--default_specie`, so it needs one. See
[Species](guides/species.md).

### `The following invalid input values have been detected`

nf-schema found a parameter of the wrong type or value. Each line names the parameter, the value
it got in brackets, and what it expected; [Parameters](reference/parameters.md) lists them all.
A `--starsolo ture` stops here rather than running with the flag off.

### Nextflow refuses the version

The manifest requires Nextflow 26.04.6 or later. On the Sanger farm the default on `PATH` is
older: `module load cellgen/nextflow/26.04.6`.

## Warnings at launch

### `Module cellgeni/... has local modifications (checksum mismatch)`

A registry module under `modules/cellgeni/` differs from the version published to the Nextflow
module registry, as recorded in its `.module-info`. It is harmless: the local copy is what runs.
It clears once the change is published (see [Developers](developers.md#modules)).

## Tasks are refused by the executor

### `Process requirement exceeds available memory` (or `CPUs`)

The local executor will not start a task bigger than the machine: the heavy steps ask for 16
cpus, and STARsolo for 48 GB or more. Cap the requests with `process.resourceLimits` in a `-c`
config — see [Resource caps](guides/profiles.md#resource-caps). Below about 32 GB of RAM a human
STAR index cannot be loaded at all.

### `Cannot exceed queue's hard limit(s). Job not submitted.`

LSF refused a request bigger than the queue allows. Under `-profile sanger`, requests above
683 GB already go to `hugemem`; if you see this, a `-c` config has set a queue by hand. On another
cluster, make your queue rule send big requests to a queue that takes them.

### Docker on Apple silicon

The images are built for `linux/amd64`. Add `--platform linux/amd64` to `docker.runOptions` in a
`-c` config; expect emulation to be slow.

## Samples go missing

### `Sample ... has no usable species '...' — skipping alignment`

The metadata named no species the pipeline can align, and no `--default_specie` was set. The
sample is listed in `results/skipped_samples.tsv`. Set `--default_specie` if you know the species
— see [Species](guides/species.md).

### `Sample ... has unknown species — using 'Homo sapiens'`

Not an error: the metadata had no species and `--default_specie` filled it in. Check the sample
if the default might be wrong for it.

### `Sample ... has runs with differing chemistries [...] — leaving chemistry to the aligners' own detection`

The sample's runs were inferred as different chemistries, so none is forced; each aligner
detects its own. This happens, for example, when a sample's runs were sequenced with different
kits, or include a library that is not gene expression. See
[Chemistry inference](guides/chemistry.md).

### RENAME10XRUN fails with `No supported 10x run layout matched whitelist evidence plus read geometry`

The run's barcode read matched no whitelist in a geometry the pipeline can align. Most often the
data is not single-cell gene expression — scATAC, a feature-barcode or V(D)J library, or not 10x
at all — and the rejection is right. Check [Archive pathologies](archive-pathologies.md) and the
[10x chemistry reference](10x_chemistry_reference.md) before concluding the submission is wrong.

### `No known 10x whitelist files found in ...`

`--wl_basedir` points at a directory without the 10x whitelists under their Cell Ranger names.
Point it at the copied `whitelists/` directory itself, not its parent.

### The run finished, but samples have no output

Failures are ignored once a task has used its retries, so the run's exit status says nothing
about lost samples. Open the run report's **Samples without a mapped output** table
([Run report](reference/run-report.md)); with `-profile lineage_keep`, `cleanup/rerun.tsv` is a
ready-made `--datasets` table of them ([Reruns and cleanup](guides/reruns.md)).

## Disk

### `Disk quota exceeded`

Raw data fills disks fast. Free space, or use `-profile lineage_keep` with `--no_raw` so
finished samples' FASTQs go as the run proceeds ([Reruns and cleanup](guides/reruns.md)).

{: .warning }
A STARsolo task that ran while the disk was full can finish "successfully" with an empty or
partial output, and `-resume` then reuses it. After a quota outage, rerun the affected samples in
a fresh work directory rather than resuming over them.
