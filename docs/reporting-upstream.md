# Reporting a defect upstream

When reprocessing turns up something genuinely wrong in a public submission, the finding is
worth more to the archive than to us. This page is how to send it, and how to be sure it is
worth sending.

It covers **archives** — GEO, SRA, ENA, ArrayExpress. These are helpdesks, not repositories:
there is no issue tracker, no maintainer reading a diff, and the evidence that convinces them is
not the evidence that convinces a developer. For bugs in this pipeline, open an issue on
[cellgeni/nf-reprocessing-public-10x](https://github.com/cellgeni/nf-reprocessing-public-10x/issues)
instead.

## 1. Check it is actually a defect

Most of what a batch rejects is not. Before anything else, read the row's table in
[Archive pathologies](archive-pathologies.md) and answer one question:

> **Would the depositor or the archive have to change something?**

If the submission is accurate and we simply downloaded it without reading `library_strategy`,
the answer is no — that is our gap, not theirs. Roughly 597 GB of one run's rejections were
correctly labelled ATAC-seq, bisulfite-seq and bulk RNA-seq. A ticket about any of them would be
wrong on its face and would spend credibility needed for the real ones.

What *is* worth reporting:

| Finding | Why it qualifies |
|---|---|
| Two libraries sharing one BioSample, so they are indistinguishable by accession | Internally inconsistent metadata; silently corrupts anyone's reprocessing |
| Reads stored at zero length, discarded on `fastq-dump` | The archived object is unusable |
| A sample with no experiment or run ID in the SRA export | The record points at nothing |
| SDRF naming FASTQ files that do not resolve | Broken references in submitted metadata |
| Two archive sources disagreeing about which files exist | One of them is wrong |

What is not: a chemistry we cannot infer, a layout the pipeline declines to guess at, a dataset
that is simply not 10x, or anything not reproduced outside our own run.

## 2. Gather evidence in the archive's terms

A maintainer wants a stack trace; an archive wants to locate one record and see that it is
broken. Five things, and the more specific the accession the better — a run (`SRR...`) beats a
sample (`GSM...`) beats a series (`GSE...`):

1. **The accession that reproduces it**, at the most specific level available.
2. **The literal field and value**, quoted from the archive's own export — `library_strategy`,
   the BioSample id, the SDRF URI. Say which file it came from.
3. **The verbatim error**, in a fenced block, from the tool that hit it. Not paraphrased: the
   string is how the helpdesk finds the record and how the next person finds the ticket.
4. **How it was measured**, with the command. "651 distinct 30-mers in the first 50,000 reads of
   R2, top one at 82.8%" is checkable; "looks like a tag library" is not.
5. **The scope** — how many runs, samples or bytes are affected. This is what moves a ticket up
   a queue, and usually the only part that makes a case for urgency.

Skip a suggested fix unless you are sure of it, and skip it *loudly* — "we have not checked
whether the other 22 samples share the BioSample" — rather than guessing and being wrong.

### The shape

````markdown
**Accession.** GSE247111 (GEO), runs SRR26669510 and SRR26669511.

**What is wrong.** Each gene-expression library and its CellPlex multiplexing-tag
library are registered against the same BioSample, so the two GSMs resolve to an
identical SRS. Any pipeline mapping runs to samples by BioSample merges a tag
library into a gene-expression sample without an error.

**Evidence.** R2 sequence complexity over the first 50,000 reads of each run:

```
SRR26669510    651 distinct 30-mers, most common 82.8%   -> tag library
SRR26669511 15,170 distinct 30-mers, most common  5.9%   -> gene expression
```

Downstream, merged outputs show 32-39% valid barcodes against 97.9% for an
unaffected sample processed in the same run.

**Scope.** Affects all 22 samples of the series; 10 of our outputs had to be
withdrawn.

**What we have not checked.** Whether the pairing is deliberate on the
depositor's side; we confirmed the shared BioSample directly for one pair and
inferred it for the other nine from an identical signature.
````

## 3. Where to send it

> **Confirm the current route on the archive's own contact page before sending.** Helpdesk
> addresses change and are not reliably discoverable from documentation mirrors. The entry
> points below are stable; read the specific address off the page when you file. Automated
> checks of the NCBI and EBI contact pages were blocked when this page was written, so **none of
> the addresses named here has been verified against a primary source.**

| Archive | Start here | Notes |
|---|---|---|
| GEO | [GEO info pages](https://www.ncbi.nlm.nih.gov/geo/info/) | GEO curators handle series and sample records. `geo@ncbi.nlm.nih.gov` is the long-standing curator address and `info@ncbi.nlm.nih.gov` is NCBI general help — confirm which is current |
| SRA | [NCBI Support Center](https://support.nlm.nih.gov/) | For run-level and archived-object problems: zero-length reads, dumps that do not match the declared layout |
| ENA | [ENA Browser support form](https://www.ebi.ac.uk/ena/browser/support) | The documented route is the support form in the ENA Browser, not an address |
| ArrayExpress / BioStudies | [BioStudies](https://www.ebi.ac.uk/biostudies/), then the EMBL-EBI contact route | SDRF and file-listing problems. Say explicitly which of the two sources you read |

A GEO record and its SRA runs are curated separately. A metadata problem in the series goes to
GEO; a broken archived object goes to SRA. Sending either to the wrong desk costs a round trip.

## 4. Draft, show, then send

Filing is outward-facing, carries the group's name, and is the one step here that cannot be
taken back quietly.

- **Draft it next to the evidence first**, in the repo or the ticket folder, so error text gets
  copied rather than retyped.
- **Check the record for an existing erratum or announcement** — a known, already-corrected
  problem costs the helpdesk more than silence.
- **An agent drafts and shows; a person sends.** This is not optional. No automated step should
  put a message into an archive's helpdesk queue.

Once filed, put the reference in the `Reported` column of
[Archive pathologies](archive-pathologies.md), next to the row it came from. That closes the
loop: the next batch that hits the same accession finds the report instead of rediagnosing it.
