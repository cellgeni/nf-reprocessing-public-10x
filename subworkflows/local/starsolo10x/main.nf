include { STARSOLO10X as REPROCESS10X_STARSOLO10X } from 'cellgeni/starsolo10x'
workflow STARSOLO10X {

    take:
    fastq_files   // channel: [ val(meta), [ file(fastq) ] ] meta: [ id: sample_id, dataset_id: dataset_id ]
    reference     // channel: [ val(meta), path(genome_reference) ] meta: [ id: reference_id ]

    main:

    // Predict memory, run limit and queue from the sample's estimated read count
    // (lib/StarsoloResources.groovy) and carry them in meta.resources, which
    // configs/starsolo10x.config reads. The estimate has to happen here: a config
    // directive closure cannot see lib/ classes in Nextflow 26. A null plan (the
    // estimate failed) makes the config fall back to the fixed ladder on `long`.
    with_resources = fastq_files
        .combine(reference.map { ref_meta, _ref -> ref_meta.id })
        .map { meta, fastqs, ref_id ->
            tuple(meta + [resources: StarsoloResources.plan(fastqs, ref_id)], fastqs)
        }

    // Run STARsolo. QC is collected in workflow/main.nf rather than here: this
    // subworkflow is instantiated once per species, whereas MAPPINGQC runs once
    // per dataset across every species it contains.
    REPROCESS10X_STARSOLO10X(with_resources, reference)

    emit:
    mapping      = REPROCESS10X_STARSOLO10X.out.mapping
    versions     = REPROCESS10X_STARSOLO10X.out.versions.first()

}
