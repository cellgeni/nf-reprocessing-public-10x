include { REPROCESS10X_STARSOLO } from '../../../modules/local/reprocess10x/starsolo'
include { REPROCESS10X_MAPPINGQC } from '../../../modules/local/reprocess10x/mappingqc'
workflow STARSOLO10X {

    take:
    fastq_files   // channel: [ val(meta), [ file(fastq) ] ] meta: [ id: sample_id, dataset_id: dataset_id ]
    reference     // channel: path(genome_reference)
    wl_basedir    // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)

    main:

    // Run STARsolo
    REPROCESS10X_STARSOLO(fastq_files, reference, wl_basedir)

    // Group by dataset
    samples_by_dataset = REPROCESS10X_STARSOLO.out.mapping
        .view { it -> "MAPPING: $it" }
        .map { meta, sample_dir -> [ groupKey([id: meta.dataset_id.getGroupTarget()], meta.dataset_id.getGroupSize()), sample_dir ] }
        .view { it -> "SAMPLES BY DATASET: $it" }
        .groupTuple()
        .view { it -> "SAMPLES BY DATASET GROUPED: $it" }
    
    // Collect mapping QC stats
    REPROCESS10X_MAPPINGQC(samples_by_dataset)

    // Collect versions
    versions = REPROCESS10X_STARSOLO.out.versions.first()
        .mix(REPROCESS10X_MAPPINGQC.out.versions.first())

    emit:
    mapping      = REPROCESS10X_STARSOLO.out.mapping
    qc_stats     = REPROCESS10X_MAPPINGQC.out.tsv
    versions     = versions

}