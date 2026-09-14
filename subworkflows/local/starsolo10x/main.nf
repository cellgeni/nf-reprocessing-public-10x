include { STARSOLO10X as REPROCESS10X_STARSOLO10X } from 'cellgeni/starsolo10x'
workflow STARSOLO10X {

    take:
    fastq_files   // channel: [ val(meta), [ file(fastq) ] ] meta: [ id: sample_id, dataset_id: dataset_id ]
    reference     // channel: [ val(meta), path(genome_reference) ] meta: [ id: reference_id ]

    main:

    // Run STARsolo. QC is collected in workflow/main.nf rather than here: this
    // subworkflow is instantiated once per species, whereas MAPPINGQC runs once
    // per dataset across every species it contains.
    REPROCESS10X_STARSOLO10X(fastq_files, reference)

    emit:
    mapping      = REPROCESS10X_STARSOLO10X.out.mapping
    versions     = REPROCESS10X_STARSOLO10X.out.versions.first()

}
