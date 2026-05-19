// IMPORT SUBWORKFLOW
include { DOWNLOAD10X } from '../subworkflows/local/download10x/'
include { FETCH10XMETA } from 'cellgeni/fetch10xmeta'
include { STARSOLO10X as STARSOLO10X_HUMAN } from '../subworkflows/local/starsolo10x/'
include { STARSOLO10X as STARSOLO10X_MOUSE } from '../subworkflows/local/starsolo10x/'

workflow REPROCESS10X {
    take:
    datasetlist      // channel: [ val(meta), [ sample_ids ] ]
    wl_basedir       // channel: [ dirpath ] a path to whitelist base directory
    human_reference  // channel: [ tuple( [id: "human"], file(human_reference) ) ]
    mouse_reference  // channel: [ tuple( [id: "mouse"], file(mouse_reference) ) ]
    metaonlyflag     // channel: [ val(metaonlyflag) ] a flag to indicate whether to only fetch metadata without downloading data or running STARsolo (e.g. for testing or debugging)
    no_infer_specie  // channel: [ val(no_infer_specie) ] a flag to indicate whether to infer specie from metadata or not; if true, all samples will be assigned the default_specie (or 'UNKNOWN' if default_specie is not set)
    default_specie   // channel: [ val(default_specie) ] the default specie to assign to samples with unknown species in metadata; if not set, 'UNKNOWN' will be used as default specie
    rawonlyflag      // channel: [ val(rawonlyflag) ] a flag to indicate whether to only download raw data (fastq files) without running STARsolo; if true, STARsolo steps will be skipped

    main:
    // STEP 0.1: Init channels
    bams            = channel.empty()
    sras            = channel.empty()
    versions        = channel.empty()
    metadata        = channel.empty()
    resolved_fastqs = channel.empty()
    starsolo        = channel.empty()
    soloqc          = channel.empty()

    // STEP 0.2: Convert dataset list to channel
    datasets = datasetlist
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_list =  row.sample_id.split(',')
            [
                [id: groupKey(row.dataset_id, sample_list.size())],
                row.sample_id
            ]
        }
    // Step 1: Fetch metadata
    FETCH10XMETA(datasets)

    // Collect metadata files per dataset
    metadata = metadata
        .mix(
            FETCH10XMETA.out.links,
            FETCH10XMETA.out.list,
            FETCH10XMETA.out.tsv,
            FETCH10XMETA.out.txt,
            FETCH10XMETA.out.soft
        )
        .groupTuple(sort: 'hash')
        .map { meta, files -> tuple( meta, files.flatten() ) }
    
    versions = versions.mix(FETCH10XMETA.out.versions)

    // STEP 2: Download datasets
    if (!metaonlyflag) {
        DOWNLOAD10X(
            FETCH10XMETA.out.links,
            wl_basedir
        )

        // Collect fastq files per sample
        def unknown_values = [null, '', 'NULL', 'UNKNOWN']
        def specie_map = [
            human: 'Homo sapiens',
            mouse: 'Mus musculus',
        ]
        resolved_fastqs = DOWNLOAD10X.out.fastq
            .map { meta, fastqs ->
                def effective_specie
                if (no_infer_specie) {
                    effective_specie = default_specie ? specie_map.get(default_specie) : 'UNKNOWN'
                } else {
                    if (meta.specie in unknown_values) {
                        effective_specie = default_specie ? specie_map.get(default_specie) : 'UNKNOWN'
                        log.warn "Sample ${meta.id} (${meta.dataset_id}) has unknown species — using effective_specie='${effective_specie}'"
                    } else if (meta.specie in specie_map.values()) {
                        effective_specie = meta.specie
                    } else {
                        effective_specie = 'UNKNOWN'
                    }
                }
                [meta + [specie: effective_specie], fastqs]
            }

        // Group samples by specie
        fastqs = resolved_fastqs
            .branch { meta, _fastqs ->
                human: meta.specie == 'Homo sapiens'
                mouse: meta.specie == 'Mus musculus'
                other: true
            }

        fastqs.other
            .map { meta, _fastqs ->
                log.warn "Sample ${meta.id} (${meta.dataset_id}) has unexpected species '${meta.specie}' — skipping STARsolo"
                [meta, _fastqs]
            }
        
        // Collect channels
        bams     = bams.mix(DOWNLOAD10X.out.bam)
        sras     = sras.mix(DOWNLOAD10X.out.sra)
        versions = versions.mix(DOWNLOAD10X.out.versions)
    }

    // STEP 3: Run STARsolo
    if (!metaonlyflag && !rawonlyflag) {
        // Run STARsolo on fastq files for human and mouse samples
        STARSOLO10X_HUMAN(fastqs.human, human_reference, wl_basedir)
        STARSOLO10X_MOUSE(fastqs.mouse, mouse_reference, wl_basedir)

        // Collect outputs
        starsolo = starsolo.mix(
            STARSOLO10X_HUMAN.out.mapping,
            STARSOLO10X_MOUSE.out.mapping
        )
        soloqc   = soloqc.mix(
            STARSOLO10X_HUMAN.out.qc_stats,
            STARSOLO10X_MOUSE.out.qc_stats
        )
        versions = versions
            .mix(
                STARSOLO10X_HUMAN.out.versions,
                STARSOLO10X_MOUSE.out.versions
            )
    }
    
    emit:
    metadata = metadata
    fastq    = resolved_fastqs
    bam      = bams
    sra      = sras
    starsolo = starsolo
    soloqc   = soloqc
    versions = versions
}