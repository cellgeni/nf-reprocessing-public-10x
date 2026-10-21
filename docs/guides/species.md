---
title: Species
parent: Guides
nav_order: 3
---

# Species

Each sample resolves to `Homo sapiens`, `Mus musculus` or `UNKNOWN`. Only the first two reach an
aligner, each against its own reference; `UNKNOWN` samples are logged, listed in
`results/skipped_samples.tsv`, and dropped.

By default the species is read from the archive metadata. A sample's runs often carry missing or
conflicting annotations, so the whole of a series' metadata is read before any sample is
decided. When a sample's species is blank, unrecognised or contradictory:

- with `--default_specie human` (or `mouse`), that species is used and a warning names the
  sample;
- without it, the sample is `UNKNOWN` and is skipped.

`--no_infer_specie` ignores the metadata altogether and gives every sample `--default_specie`,
which it then requires. Use it when you know the whole table is one species and the metadata is
not to be trusted.

{: .note }
A species is aligned only if its reference is set: `--star_human_reference` /
`--star_mouse_reference` for STARsolo, `--cr_human_reference` / `--cr_mouse_reference` for Cell
Ranger. `-profile sanger` sets all four; elsewhere, a species without one is not aligned and a
warning at launch says so.
