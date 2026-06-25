include { REPROCESS10X_LOADDATA } from '../../../modules/local/reprocess10x/loaddata'
include { REPROCESS10X_BAM2FASTQ } from '../../../modules/local/reprocess10x/bam2fastq'
include { REPROCESS10X_SRA2FASTQ } from '../../../modules/local/reprocess10x/sra2fastq'
include { REPROCESS10X_RENAMEFASTQ } from '../../../modules/local/reprocess10x/renamefastq'

workflow DOWNLOAD10X {

    take:
    links     // channel: [ val(meta), path("links.tsv") ]
    wl_basedir   // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)
    
    main:
    // STEP 0: Initialize channels
    sras     = channel.empty()
    bams     = channel.empty()
    versions = channel.empty()

    // STEP 1: Load data from links
    // Get links for each run file from metadata
    collected_links = links
        // Read links file and split by tab
        .splitCsv(sep: '\t', strip: true)
        // Group by sample to count number of runs per sample
        .map { meta, row -> tuple( row[4], [meta, row] ) }
        .groupTuple(sort: 'hash', remainder: true)
        .map { _sample_id, metas_rows -> tuple( metas_rows.size(), metas_rows ) }
        .transpose()
        // Flatten the rows
        .flatMap { size, meta_row -> 
            def (meta, row) = meta_row
            def run_meta = [
                id       : row[0], // run ID
                sample_id: groupKey(row[4], size), // sample ID
                dataset_id : meta.id, // dataset's series ID
                specie   : row[1], // specie
                type     : row[3] // file type i.e. BAM, FASTQ, SRA ...
            ]
            def urls = row[2].split(";")
            return urls.collect { url -> [groupKey(run_meta, urls.size()), url] }
        }
        //.view { meta, url -> "LINKS: meta=[${meta.getGroupTarget().collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], url=$url (${url.getClass().simpleName}), count=${meta.getGroupSize()}" }
                                           

    REPROCESS10X_LOADDATA(collected_links)

    // STEP 2: Convert loaded data to fastq if needed
    bams = bams.mix(REPROCESS10X_LOADDATA.out.bam)
    REPROCESS10X_BAM2FASTQ(bams)
    
    sras = sras.mix(REPROCESS10X_LOADDATA.out.sra)
    REPROCESS10X_SRA2FASTQ(sras, wl_basedir)

    // STEP 3 Collect all outputs
    // Combine all fastq channels and group by sample
    fastqs = REPROCESS10X_LOADDATA.out.fastq
        // Combine fastq files for each read as they were loaded separately
        .groupTuple(sort: 'hash', remainder: true)
        // Combine fastq files from BAM and SRA conversion
        .mix(
            REPROCESS10X_BAM2FASTQ.out.fastq,
            REPROCESS10X_SRA2FASTQ.out.fastq
        )
        // Leave only sample id and dataset id in metadata
        .map { run_meta, fastq ->
            def sample_meta = [id: run_meta.sample_id.getGroupTarget(), dataset_id: run_meta.dataset_id, specie: run_meta.specie]
            def run_count = run_meta.sample_id.getGroupSize()
            tuple( groupKey(sample_meta, run_count), fastq )
        }
        // Group by sample id and dataset id
        .groupTuple(sort: 'hash', remainder: true)
        // Combine a list of fastq files
        .map { groupkey, fastqlist -> tuple( groupkey.getGroupTarget(), fastqlist.flatten() ) }

    REPROCESS10X_RENAMEFASTQ(fastqs, wl_basedir)

    // Collect versions
    versions = versions
        .mix(
            REPROCESS10X_LOADDATA.out.versions.first(),
            REPROCESS10X_BAM2FASTQ.out.versions.first(),
            REPROCESS10X_SRA2FASTQ.out.versions.first(),
            REPROCESS10X_RENAMEFASTQ.out.versions.first()
        )

    emit:
    fastq    = REPROCESS10X_RENAMEFASTQ.out.fastq
    fastq_inference_report  = REPROCESS10X_RENAMEFASTQ.out.report
    fastq_inference_summary = REPROCESS10X_RENAMEFASTQ.out.summary
    bam      = bams
    sra      = sras
    versions = versions
}
