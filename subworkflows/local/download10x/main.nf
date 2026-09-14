include { WGET10X } from '../../../modules/cellgeni/wget10x'
include { REPROCESS10X_BAM2FASTQ } from '../../../modules/local/reprocess10x/bam2fastq'
include { SRA2FASTQ } from '../../../modules/cellgeni/sra2fastq'
include { RENAME10XRUN } from '../../../modules/cellgeni/rename10xrun'
include { RENAME10XSAMPLE } from '../../../modules/cellgeni/rename10xsample'

// Settle one species for a sample from its per-run annotations, which are
// frequently inconsistent (SRA reports 'UNKNOWN' for a subset of the runs, and
// occasionally two different organisms for the same sample).
//
//   single known species  -> that species; the sample's UNKNOWN runs inherit it
//   several known species -> unresolvable, falls back to --default_specie
//   all unknown           -> --default_specie
//   a species we have no reference for -> left UNKNOWN, so the sample is
//                            skipped rather than forced onto the default genome
//
// Returns 'UNKNOWN' when nothing can be decided; the species branch in
// workflow/main.nf drops those samples.
def resolve_specie(sample_id, dataset_id, run_species, no_infer_specie, default_specie) {
    def specie_map     = [human: 'Homo sapiens', mouse: 'Mus musculus']
    def unknown_values = [null, '', 'NULL', 'UNKNOWN']
    def fallback       = default_specie ? specie_map.get(default_specie.toLowerCase(), 'UNKNOWN') : 'UNKNOWN'

    if (no_infer_specie) {
        return fallback
    }

    def known = run_species.unique().findAll { specie -> !(specie in unknown_values) }

    if (known.size() > 1) {
        log.warn "Sample ${sample_id} (${dataset_id}) has conflicting species across its runs ${known} — using '${fallback}'"
        return fallback
    }
    if (known.isEmpty()) {
        log.warn "Sample ${sample_id} (${dataset_id}) has unknown species — using '${fallback}'"
        return fallback
    }
    if (!(known[0] in specie_map.values())) {
        log.warn "Sample ${sample_id} (${dataset_id}) has species '${known[0]}', which has no reference"
        return 'UNKNOWN'
    }
    return known[0]
}

// Settle one chemistry for a sample from the per-run reports RENAME10XRUN wrote, or
// null when there is nothing to force. Both aligners take their chemistry from this
// single answer, so they can never be handed contradicting ones.
def resolve_run_chemistry(sample_id, chemistry_files) {
    def chemistry_ids = chemistry_files
        .findAll { report -> report }
        .collect { report ->
            try {
                new groovy.json.JsonSlurper().parseText(report.text)
            } catch (Exception err) {
                log.warn "Sample ${sample_id}: could not read chemistry report ${report.name} (${err.message})"
                return null
            }
        }
        .findAll { report -> report && report.ok && report.chemistry_id }
        .collect { report -> report.chemistry_id }
        .unique()

    // Runs that disagree are not safe to force onto one chemistry, and a sample with no
    // usable report (every BAM-derived one, for instance) has nothing to say either way.
    if (chemistry_ids.size() > 1) {
        log.warn "Sample ${sample_id} has runs with differing chemistries ${chemistry_ids} — leaving chemistry to the aligners' own detection"
        return null
    }

    return chemistry_ids.isEmpty() ? null : chemistry_ids[0]
}

// Turn an inferred chemistry id into the value Cell Ranger's --chemistry expects,
// or null to leave it on auto-detection.
//
// Cell Ranger re-detects chemistry from the FASTQs, and for Multiome gene expression it
// detects ARC-v1 and then refuses to proceed unless the chemistry was named explicitly —
// so a whole Multiome dataset fails at the aligner even though run-level inference had
// already identified it correctly. Passing what we inferred settles that.
//
// Only chemistries that are both called unambiguously by RENAME10XRUN and not already
// handled by --chemistry=auto are mapped. gex_3pv2_or_5pv1v2 and gex_3pv3_family are
// deliberately absent: infer_10x_run marks them 'layout_only' because one whitelist
// covers several chemistries, and auto resolves them correctly today. Forcing a guess
// there would risk the samples that currently work in exchange for nothing.
def resolve_cellranger_chemistry(chemistry_id) {
    def cellranger_chemistry = [ 'gex_multiome_arc_v1': 'ARC-v1' ]

    return chemistry_id ? cellranger_chemistry[chemistry_id] : null
}

// Turn an inferred chemistry id into the value `starsolo 10x --wl` expects, or null to
// leave the wrapper to detect the whitelist itself.
//
// The wrapper accepts infer_10x_run.py's own chemistry ids, so a resolved id goes
// straight through. Unlike Cell Ranger, the layout_only ids are safe to pass here: all
// --wl decides is the whitelist and the CB/UMI geometry, and the wrapper derives that
// geometry from the whitelist file alone — so a forced id gives exactly what its own
// detection would have given, minus the detection. What that avoids is the failure mode:
// the wrapper matches barcodes from a 200k subsample and aborts when no whitelist clears
// its threshold, on data that matched a whitelist perfectly well at run level.
//
// ATAC ids are rejected by the wrapper (`starsolo 10x` is gene expression only), so they
// are dropped here rather than turned into a task that dies on its first line.
def resolve_starsolo_chemistry(sample_id, chemistry_id) {
    def starsolo_chemistries = [
        'gex_3pv1',
        'gex_3pv2_or_5pv1v2',
        'gex_3pv3_family',
        'gex_3pv4_gemx',
        'gex_5pv3_gemx',
        'gex_multiome_arc_v1',
    ]

    if (!chemistry_id) {
        return null
    }
    if (!(chemistry_id in starsolo_chemistries)) {
        log.warn "Sample ${sample_id} was called as '${chemistry_id}', which `starsolo 10x --wl` does not accept — leaving STARsolo to detect the whitelist"
        return null
    }

    return chemistry_id
}

workflow DOWNLOAD10X {

    take:
    links     // channel: [ val(meta), path("links.tsv") ]
    wl_basedir   // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)
    no_infer_specie  // value: skip reading species from metadata; assign default_specie to every sample
    default_specie   // value: species to assign when metadata is missing, unknown or contradictory

    main:
    // STEP 0: Initialize channels
    sras     = channel.empty()
    bams     = channel.empty()
    versions = channel.empty()

    // STEP 1: Load data from links
    // links.tsv holds one row per run, one file per dataset. It is read whole
    // rather than row by row so that species, per-sample run counts and the
    // dataset's sample count are all settled before any of them is used as a
    // grouping key — none of these can be decided from a single run's row.
    collected_links = links
        .flatMap { meta, links_file ->
            def rows       = links_file.splitCsv(sep: '\t', strip: true)
            def dataset_id = meta.id.toString()
            def by_sample  = rows.groupBy { row -> row[4] }

            // One species per sample, decided from all of its runs at once, so
            // that every run of a sample carries the same value downstream.
            def specie_by_sample = by_sample.collectEntries { sample_id, sample_rows ->
                [ (sample_id): resolve_specie(sample_id, dataset_id, sample_rows.collect { row -> row[1] }, no_infer_specie, default_specie) ]
            }

            // Samples that will actually reach an aligner. This is what MAPPINGQC
            // waits for, so it must exclude the ones skipped for want of a species.
            def aligned_samples = specie_by_sample.count { _sample_id, specie -> specie in ['Homo sapiens', 'Mus musculus'] }

            rows.collectMany { row ->
                def run_meta = [
                    id        : row[0], // run ID
                    sample_id : groupKey(row[4], by_sample[row[4]].size()), // sample ID + its own run count
                    dataset_id: groupKey(dataset_id, aligned_samples), // series ID + its alignable sample count
                    specie    : specie_by_sample[row[4]],
                    type      : row[3] // file type i.e. BAM, FASTQ, SRA ...
                ]
                def urls = row[2].split(";")
                return urls.collect { url -> [groupKey(run_meta, urls.size()), url] }
            }
        }
        //.view { meta, url -> "LINKS: meta=[${meta.getGroupTarget().collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], url=$url (${url.getClass().simpleName}), count=${meta.getGroupSize()}" }
                                          

    WGET10X(collected_links)

    //REPROCESS10X_LOADDATA.out.fastq.view { meta, fastq -> "FASTQ: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastq=$fastq (${fastq.getClass().simpleName})" }
    //REPROCESS10X_LOADDATA.out.sra.view { meta, sra -> "SRA: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], sra=$sra (${sra.getClass().simpleName})" }
    //REPROCESS10X_LOADDATA.out.bam.view { meta, bam -> "BAM: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], bam=$bam (${bam.getClass().simpleName})" }

    // STEP 2: Convert loaded data to fastq if needed
    bams = bams.mix(WGET10X.out.bam)
    REPROCESS10X_BAM2FASTQ(bams)
    
    sras = sras.mix(WGET10X.out.sra)
    SRA2FASTQ(sras)

    // STEP 3: Rename fastq files to match 10x run naming convention
    fastqs2rename = WGET10X.out.fastq
        // Combine fastq files for each read as they were loaded separately
        .groupTuple(sort: 'hash', remainder: true)
        .mix(SRA2FASTQ.out.fastq)

    RENAME10XRUN(
        fastqs2rename,
        wl_basedir
    )

    // Collect every per-run fastq for each sample before renaming: renamed
    // ENA/SRA runs (RENAME10XRUN) together with BAM-derived runs
    // (REPROCESS10X_BAM2FASTQ). Mixing here means a sample whose runs are of
    // mixed origin (some BAM, some FASTQ) is grouped once and handed to a
    // single RENAME10XSAMPLE job — instead of being split across two branches.
    // Species is settled per sample in STEP 1, so it is deliberately left out of
    // the grouping key: were the runs of one sample ever to disagree on it,
    // keying on it would split them into separate groups that each produce a
    // STARsolo directory named after the sample — which then collide when
    // MAPPINGQC stages them side by side.
    // Each renamed run travels with the chemistry report that describes it, so the
    // per-sample step can settle CB/UMI geometry from what run-level inference
    // already decided rather than measuring read lengths again and disagreeing with
    // it. BAM-derived runs never pass through RENAME10XRUN and so carry no report;
    // they join with an empty list, which RENAME10XSAMPLE treats as "no metadata for
    // this run" exactly as before.
    runs_with_chemistry = RENAME10XRUN.out.reads
        .join(RENAME10XRUN.out.chemistry, failOnMismatch: true, failOnDuplicate: true)
        .mix(REPROCESS10X_BAM2FASTQ.out.fastq.map { run_meta, fastq -> tuple(run_meta, fastq, []) })

    sample2rename = runs_with_chemistry
        .map { run_meta, fastq, chemistry ->
            def sample_meta = [id: run_meta.sample_id.getGroupTarget(), dataset_id: run_meta.dataset_id]
            def run_count = run_meta.sample_id.getGroupSize()
            tuple( groupKey(sample_meta, run_count), fastq, chemistry, run_meta.specie )
        }
        .groupTuple(sort: 'hash', remainder: true)
        .map { sample_key, fastqlist, chemistrylist, species ->
            def sample_meta = sample_key.getGroupTarget() + [specie: species.first()]
            def reports     = chemistrylist.flatten()

            // Both aligner keys are added here, where the run reports are in scope, and
            // each is added only when inference actually decided something: a sample with
            // no usable report keeps the meta it always had. Note that `wl` does get set
            // for nearly every FASTQ/SRA-derived sample, which puts a new value in
            // RENAME10XSAMPLE's task hash — those tasks re-run once on the next -resume
            // even though renaming itself is unaffected by chemistry. (The aligners re-run
            // regardless: Nextflow folds a process's own source text into its task hash,
            // so editing either module invalidates it.)
            def chemistry_id = resolve_run_chemistry(sample_meta.id, reports)

            def cr_chemistry = resolve_cellranger_chemistry(chemistry_id)
            if (cr_chemistry) {
                log.info "Sample ${sample_meta.id}: passing --chemistry=${cr_chemistry} to Cell Ranger from run-level inference"
                sample_meta = sample_meta + [chemistry: cr_chemistry]
            }

            def starsolo_chemistry = resolve_starsolo_chemistry(sample_meta.id, chemistry_id)
            if (starsolo_chemistry) {
                log.info "Sample ${sample_meta.id}: passing --wl ${starsolo_chemistry} to STARsolo from run-level inference"
                sample_meta = sample_meta + [wl: starsolo_chemistry]
            }

            tuple(
                sample_meta,
                fastqlist.flatten(),
                reports
            )
        }
    RENAME10XSAMPLE(
        sample2rename
        )



    // STEP 4 Collect all outputs
    // Every sample's fastqs now come from a single RENAME10XSAMPLE job (BAM- and
    // FASTQ-origin runs were merged per sample above), so each sample is emitted
    // exactly once — no more duplicate STARsolo jobs for the same sample.
    fastqs = RENAME10XSAMPLE.out.reads

    // Collect versions
    versions = versions
        .mix(
            WGET10X.out.versions.first(),
            REPROCESS10X_BAM2FASTQ.out.versions.first(),
            RENAME10XRUN.out.versions.first(),
            RENAME10XSAMPLE.out.versions.first(),
            SRA2FASTQ.out.versions.first()
        )

    emit:
    original_fastq = fastqs2rename.mix(REPROCESS10X_BAM2FASTQ.out.fastq)
    runs           = RENAME10XRUN.out.reads.mix(REPROCESS10X_BAM2FASTQ.out.fastq)
    fastq          = fastqs
    bam            = bams
    sra            = sras
    versions       = versions
}
