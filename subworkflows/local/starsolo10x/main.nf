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
        //.view { meta, mapping_dir -> "SAMPLES: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], mapping_dir=$mapping_dir (${mapping_dir.getClass().simpleName})" }
        // Create a dataset-level grouping key
        .map { meta, sample_dir ->
            def dataset_meta = [id: meta.dataset_id.getGroupTarget()]
            def sample_count = meta.dataset_id.getGroupSize()
            tuple( groupKey(dataset_meta, sample_count), sample_dir )
        }
        //.view { meta, sample_dir -> "SAMPLES BY DATASET: meta=[${meta.hasProperty('groupTarget') ? "GroupKey(${meta.getGroupTarget()}, size=${meta.getGroupSize()})" : meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], sample_dir=$sample_dir (${sample_dir.getClass().simpleName})" }
        // Group by dataset
        .groupTuple(sort: 'hash')
        //.view { meta, samples -> "SAMPLES GROUPED BY DATASET: meta=[${meta.hasProperty('groupTarget') ? "GroupKey(${meta.getGroupTarget()}, size=${meta.getGroupSize()})" : meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], samples=$samples (${samples.getClass().simpleName})" }
    
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