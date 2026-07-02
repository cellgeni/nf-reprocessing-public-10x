// IMPORT SUBWORKFLOW
include { DOWNLOAD10X } from '../subworkflows/local/download10x/'
include { FETCH10XMETA } from 'cellgeni/fetch10xmeta'
include { STARSOLO10X as STARSOLO10X_HUMAN } from '../subworkflows/local/starsolo10x/'
include { STARSOLO10X as STARSOLO10X_MOUSE } from '../subworkflows/local/starsolo10x/'
include { CELLRANGER_COUNT as CELLRANGER_COUNT_HUMAN } from '../modules/cellgeni/cellranger/count'
include { CELLRANGER_COUNT as CELLRANGER_COUNT_MOUSE } from '../modules/cellgeni/cellranger/count'

workflow REPROCESS10X {
    take:
    datasetlist      // channel: [ val(meta), [ sample_ids ] ]
    wl_basedir       // channel: [ dirpath ] a path to whitelist base directory
    star_human_reference  // channel: [ tuple( [id: "human"], file(human_reference) ) ]
    star_mouse_reference  // channel: [ tuple( [id: "mouse"], file(mouse_reference) ) ]
    cr_human_reference  // channel: [ tuple( [id: "human"], file(human_reference) ) ]
    cr_mouse_reference // channel: [ tuple( [id: "mouse"], file(mouse_reference) ) ]
    metaonlyflag     // channel: [ val(metaonlyflag) ] only fetch metadata, skip download and alignment
    no_infer_specie  // channel: [ val(no_infer_specie) ] skip reading species from metadata; assign default_specie to all samples
    default_specie   // channel: [ val(default_specie) ] species to assign when metadata is missing or unknown
    starsoloflag     // channel: [ val(starsoloflag) ] run STARsolo alignment after downloading
    cellrangerflag   // channel: [ val(cellrangerflag) ] run Cell Ranger alignment after downloading (not yet implemented)

    main:
    // STEP 0.1: Init channels
    bams            = channel.empty()
    sras            = channel.empty()
    versions        = channel.empty()
    metadata        = channel.empty()
    resolved_fastqs = channel.empty()
    starsolo        = channel.empty()
    soloqc          = channel.empty()
    cellranger      = channel.empty()

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
        resolved_fastqs = resolved_fastqs.mix(DOWNLOAD10X.out.fastq)
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

    // STEP 3.1: Run STARsolo
    if (!metaonlyflag && starsoloflag) {
        // Run STARsolo on fastq files for human and mouse samples
        STARSOLO10X_HUMAN(fastqs.human, star_human_reference)
        STARSOLO10X_MOUSE(fastqs.mouse, star_mouse_reference)

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

    // STEP 3.2: Run Cell Ranger
    if (!metaonlyflag && cellrangerflag) {
        // Run Cell Ranger on fastq files for human and mouse samples
        CELLRANGER_COUNT_HUMAN(fastqs.human, cr_human_reference)
        CELLRANGER_COUNT_MOUSE(fastqs.mouse, cr_mouse_reference)

        // Collect outputs
        cellranger = cellranger.mix(
            CELLRANGER_COUNT_HUMAN.out.mapping,
            CELLRANGER_COUNT_MOUSE.out.mapping
        )
        
        versions = versions
            .mix(
                CELLRANGER_COUNT_HUMAN.out.versions.first(),
                CELLRANGER_COUNT_MOUSE.out.versions.first()
            )
    }


    
    emit:
    metadata       = metadata
    original_fastq = DOWNLOAD10X.out.original_fastq
    runs           = DOWNLOAD10X.out.runs
    fastq          = resolved_fastqs
    bam            = bams
    sra            = sras
    starsolo       = starsolo
    soloqc         = soloqc
    cellranger     = cellranger
    versions       = versions
}