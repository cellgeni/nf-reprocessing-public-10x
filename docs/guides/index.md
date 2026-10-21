---
title: Guides
nav_order: 3
has_children: true
---

# Guides

Each guide covers one job you might need to do with the pipeline.

| Guide | Covers |
|---|---|
| [Profiles and running elsewhere](profiles.md) | `sanger`, `singularity`, `docker`; whitelists, references and resource caps outside the farm; a cluster config |
| [Resources and queues](resources.md) | The resource labels, STARsolo's per-sample prediction, how a queue is chosen |
| [Species](species.md) | How each sample's species is decided, `--default_specie`, `--no_infer_specie` |
| [Chemistry inference](chemistry.md) | How a run's chemistry is read from its FASTQs and what each aligner is told |
| [Reruns and cleanup](reruns.md) | Retries, `-resume`, `-profile lineage_keep`, `--no_raw` |
