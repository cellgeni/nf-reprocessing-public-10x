include { WGET10X } from '../../../modules/cellgeni/wget10x'
include { REPROCESS10X_BAM2FASTQ } from '../../../modules/local/reprocess10x/bam2fastq'
include { SRA2FASTQ } from '../../../modules/cellgeni/sra2fastq'
include { RENAME10XRUN } from '../../../modules/cellgeni/rename10xrun'
include { RENAME10XSAMPLE } from '../../../modules/cellgeni/rename10xsample'

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

    sample2rename = RENAME10XRUN.out.reads
        .map { run_meta, fastq ->
            def sample_meta = [id: run_meta.sample_id.getGroupTarget(), dataset_id: run_meta.dataset_id, specie: run_meta.specie]
            def run_count = run_meta.sample_id.getGroupSize()
            tuple( groupKey(sample_meta, run_count), fastq )
        }
        .groupTuple(sort: 'hash', remainder: true)
        .map { run_meta, fastqlist -> tuple( run_meta.getGroupTarget(), fastqlist.flatten() ) }
    RENAME10XSAMPLE(
        sample2rename
        )



    // STEP 4 Collect all outputs
    // Combine all fastq channels and group by sample
    fastqs = REPROCESS10X_BAM2FASTQ.out.fastq
        //.view { meta, fastqs -> "FASTQ MIXED: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Leave only sample id and dataset id in metadata
        .map { run_meta, fastq ->
            def sample_meta = [id: run_meta.sample_id.getGroupTarget(), dataset_id: run_meta.dataset_id, specie: run_meta.specie]
            def run_count = run_meta.sample_id.getGroupSize()
            tuple( groupKey(sample_meta, run_count), fastq )
        }
        //.view { groupkey, fastqs -> "FASTQ PRE-GROUPED: groupkey=$groupkey (${groupkey.getClass().simpleName}), fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Group by sample id and dataset id
        .groupTuple(sort: 'hash', remainder: true)
        //.view { groupkey, fastqs -> "FASTQ GROUPED 2: groupkey=$groupkey (${groupkey.getClass().simpleName}), fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Combine a list of fastq files
        .map { groupkey, fastqlist -> tuple( groupkey.getGroupTarget(), fastqlist.flatten() ) }
        .mix(RENAME10XSAMPLE.out.reads)
        //.view { meta, fastqs -> "FASTQ FINAL: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastqs (${fastqs.getClass().simpleName})" }

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
