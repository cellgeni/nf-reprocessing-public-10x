// IMPORT SUBWORKFLOW
include { DOWNLOAD10X } from '../subworkflows/local/download10x/'
include { FETCH10XMETA } from 'cellgeni/fetch10xmeta'
include { STARSOLO10X as STARSOLO10X_HUMAN } from '../subworkflows/local/starsolo10x/'
include { STARSOLO10X as STARSOLO10X_MOUSE } from '../subworkflows/local/starsolo10x/'
include { CELLRANGER_COUNT as CELLRANGER_COUNT_HUMAN } from '../modules/cellgeni/cellranger/count'
include { CELLRANGER_COUNT as CELLRANGER_COUNT_MOUSE } from '../modules/cellgeni/cellranger/count'
include { REPROCESS10X_MAPPINGQC } from '../modules/local/reprocess10x/mappingqc'

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
    original_fastq  = channel.empty()
    runs            = channel.empty()
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
            def sample_list =  row.containsKey("sample_id") && row.sample_id ? row.sample_id.split(',') : []
            [
                [id: sample_list.size() > 0 ? groupKey(row.dataset_id, sample_list.size()) : row.dataset_id],
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
            wl_basedir,
            no_infer_specie,
            default_specie
        )

        // Species is already resolved per sample inside DOWNLOAD10X, where the
        // whole links.tsv is in scope — a sample's runs often carry conflicting
        // or missing annotations that cannot be settled one run at a time.
        resolved_fastqs = resolved_fastqs.mix(DOWNLOAD10X.out.fastq)

        // Group samples by specie
        fastqs = resolved_fastqs
            .branch { meta, _fastqs ->
                human: meta.specie == 'Homo sapiens'
                mouse: meta.specie == 'Mus musculus'
                other: true
            }

        fastqs.other
            .map { meta, _fastqs ->
                log.warn "Sample ${meta.id} (${meta.dataset_id}) has no usable species — skipping alignment"
                [meta, _fastqs]
            }
        
        // Collect channels
        original_fastq = original_fastq.mix(DOWNLOAD10X.out.original_fastq)
        runs           = runs.mix(DOWNLOAD10X.out.runs)
        bams           = bams.mix(DOWNLOAD10X.out.bam)
        sras           = sras.mix(DOWNLOAD10X.out.sra)
        versions       = versions.mix(DOWNLOAD10X.out.versions)
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

        // Collect mapping QC stats once per dataset, over every species it
        // holds: both aligner branches carry the same dataset id, so QC'ing them
        // separately would have each publish its own <dataset>.solo_qc.tsv to
        // the same path. dataset_id's group size is the dataset's alignable
        // sample count, set in DOWNLOAD10X once species are resolved.
        samples_by_dataset = starsolo
            .map { meta, sample_dir ->
                def dataset_meta = [id: meta.dataset_id.getGroupTarget()]
                tuple( groupKey(dataset_meta, meta.dataset_id.getGroupSize()), sample_dir )
            }
            .groupTuple(sort: 'hash', remainder: true)

        REPROCESS10X_MAPPINGQC(samples_by_dataset)

        soloqc   = soloqc.mix(REPROCESS10X_MAPPINGQC.out.tsv)
        versions = versions
            .mix(
                STARSOLO10X_HUMAN.out.versions,
                STARSOLO10X_MOUSE.out.versions,
                REPROCESS10X_MAPPINGQC.out.versions.first()
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
    original_fastq = original_fastq
    runs           = runs
    fastq          = resolved_fastqs
    bam            = bams
    sra            = sras
    starsolo       = starsolo
    soloqc         = soloqc
    cellranger     = cellranger
    versions       = versions
}