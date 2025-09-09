include { REPROCESS10X_PARSEMETADATA } from '../../../modules/local/reprocess10x/parsemetadata'
include { REPROCESS10X_LOADDATA } from '../../../modules/local/reprocess10x/loaddata'
include { REPROCESS10X_BAM2FASTQ } from '../../../modules/local/reprocess10x/bam2fastq'
include { REPROCESS10X_SRA2FASTQ } from '../../../modules/local/reprocess10x/sra2fastq'

workflow DOWNLOAD10X {

    take:
    datasets     // channel: [ val(meta), [ sample_ids ] ]
    wl_basedir   // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)
    
    main:
    // STEP 1: Parse metadata for each dataset
    REPROCESS10X_PARSEMETADATA(datasets)

    // Get links for each run file from metadata
    links = REPROCESS10X_PARSEMETADATA.out.links
        // Read links file and split by tab
        .splitCsv(sep: '\t', strip: true)
        // Group by sample to count number of runs per sample
        .map { meta, row -> tuple( row[4], [meta, row] ) }
        .groupTuple(sort: 'hash')
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
            return urls.collect { url -> [run_meta, url] }
        }
        //.view { meta, url -> "LINKS: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], url=$url (${url.getClass().simpleName})" }
                                          

      
    // STEP 2: Load data from links
    REPROCESS10X_LOADDATA(links)

    //REPROCESS10X_LOADDATA.out.fastq.view { meta, fastq -> "FASTQ: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastq=$fastq (${fastq.getClass().simpleName})" }
    //REPROCESS10X_LOADDATA.out.sra.view { meta, sra -> "SRA: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], sra=$sra (${sra.getClass().simpleName})" }
    //REPROCESS10X_LOADDATA.out.bam.view { meta, bam -> "BAM: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], bam=$bam (${bam.getClass().simpleName})" }

    // // STEP 3: Convert loaded data to fastq if needed
    REPROCESS10X_BAM2FASTQ(REPROCESS10X_LOADDATA.out.bam)
    REPROCESS10X_SRA2FASTQ(REPROCESS10X_LOADDATA.out.sra, wl_basedir)

    // STEP 4: Collect all outputs
    // Combine all fastq channels and group by sample
    fastqs = REPROCESS10X_LOADDATA.out.fastq
        // Combine fastq files for each run as they were loaded separately
        //.view { meta, fastq -> "FASTQ LOADED: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastq=$fastq (${fastq.getClass().simpleName})" }
        .groupTuple(size: 2, sort: 'hash')
        //.view { meta, fastqs -> "FASTQ GROUPED: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Combine fastq files from BAM and SRA conversion
        .mix(
            REPROCESS10X_BAM2FASTQ.out.fastq,
            REPROCESS10X_SRA2FASTQ.out.fastq
        )
        //.view { meta, fastqs -> "FASTQ MIXED: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Leave only sample id and dataset id in metadata
        .map { run_meta, fastq ->
            def sample_meta = [id: run_meta.sample_id.getGroupTarget(), dataset_id: run_meta.dataset_id, specie: run_meta.specie]
            def run_count = run_meta.sample_id.getGroupSize()
            tuple( groupKey(sample_meta, run_count), fastq )
        }
        //.view { groupkey, meta, fastqs -> "FASTQ PRE-GROUPED: groupkey=$groupkey (${groupkey.getClass().simpleName}), meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Group by sample id and dataset id
        .groupTuple(sort: 'hash')
        //.view { groupkey, fastqs -> "FASTQ GROUPED 2: groupkey=$groupkey (${groupkey.getClass().simpleName}), fastqs=$fastqs (${fastqs.getClass().simpleName})" }
        // Combine a list of fastq files
        .map { groupkey, fastqlist -> tuple( groupkey.getGroupTarget(), fastqlist.flatten() ) }
        //.view { meta, fastqs -> "FASTQ FINAL: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastqs (${fastqs.getClass().simpleName})" }
    
    // Collect metadata files per dataset
    metadata = REPROCESS10X_PARSEMETADATA.out.links
        .mix(
            REPROCESS10X_PARSEMETADATA.out.list,
            REPROCESS10X_PARSEMETADATA.out.tsv,
            REPROCESS10X_PARSEMETADATA.out.txt,
            REPROCESS10X_PARSEMETADATA.out.soft
        )
        .groupTuple(sort: 'hash')
        .map { meta, files -> tuple( meta, files.flatten() ) }
    
    // Collect versions
    versions = REPROCESS10X_PARSEMETADATA.out.versions.first()
        .mix(
            REPROCESS10X_LOADDATA.out.versions.first(),
            REPROCESS10X_BAM2FASTQ.out.versions.first(),
            REPROCESS10X_SRA2FASTQ.out.versions.first()
        )

    emit:
    fastq    = fastqs
    metadata = metadata
    versions = versions
}
