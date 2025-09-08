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
        .groupTuple()
        .map { _sample_id, metas_rows -> tuple( metas_rows.size(), metas_rows ) }
        .transpose()
        // Flatten the rows
        .flatMap { size, meta_row -> 
            def (meta, row) = meta_row
            def run_meta = [
                id       : row[0], // run ID
                sample_id: row[4], // sample ID
                sample_run_counts: size, // number of runs for this sample
                dataset_id : meta.id, // dataset's series ID
                dataset_sample_count: meta.dataset_sample_count, // number of samples in this dataset
                specie   : row[1], // specie
                type     : row[3] // file type i.e. BAM, FASTQ, SRA ...
            ]
            def urls = row[2].split(";")
            return urls.collect { url -> [run_meta, url] }
        }
        .view { it -> "LINKS: $it" }
                                          

      
    // STEP 2: Load data from links
    REPROCESS10X_LOADDATA(links)

    REPROCESS10X_LOADDATA.out.fastq.view { it -> "FASTQ: $it" }
    REPROCESS10X_LOADDATA.out.sra.view() { it -> "SRA: $it" }
    REPROCESS10X_LOADDATA.out.bam.view() { it -> "BAM: $it" }

    // STEP 3: Convert loaded data to fastq if needed
    REPROCESS10X_BAM2FASTQ(REPROCESS10X_LOADDATA.out.bam)
    REPROCESS10X_SRA2FASTQ(REPROCESS10X_LOADDATA.out.sra, wl_basedir)

    // STEP 4: Collect all outputs
    // Combine all fastq channels and group by sample
    fastq_files = REPROCESS10X_LOADDATA.out.fastq
        // Combine fastq files for each run as they were loaded separately
        .view { it -> "FASTQ LOADED: $it" }
        .groupTuple(size: 2)
        .view { it -> "FASTQ GROUPED: $it" }
        // Combine fastq files from BAM and SRA conversion
        .mix(
            REPROCESS10X_BAM2FASTQ.out.fastq,
            REPROCESS10X_SRA2FASTQ.out.fastq
        )
        .view { it -> "FASTQ MIXED: $it" }
        // Leave only sample id and dataset id in metadata
        .map { run_meta, fastq ->
            [
                run_meta.sample_id,
                [
                    id       : run_meta.sample_id.getGroupTarget(),
                    dataset_id : run_meta.dataset_id,
                    specie: run_meta.specie,
                    type: run_meta.type
                ],
                fastq
            ]
        }
        .view { it -> "FASTQ PRE-GROUP: $it" }
        // Group by sample id and dataset id
        .groupTuple()
        .view { it -> "FASTQ GROUPED 2: $it" }
        // Combine a list of fastq files
        .map { _groupkey, meta, fastqlists -> tuple( meta, fastqlists.flatten() ) }
        .view { it -> "FASTQ FINAL: $it" }
    
    // Collect metadata files per dataset
    metadata = REPROCESS10X_PARSEMETADATA.out.links
        .mix(
            REPROCESS10X_PARSEMETADATA.out.list,
            REPROCESS10X_PARSEMETADATA.out.tsv,
            REPROCESS10X_PARSEMETADATA.out.txt,
            REPROCESS10X_PARSEMETADATA.out.soft
        )
        .groupTuple()
        .map { meta, files -> tuple( meta, files.flatten() ) }
    
    // Collect versions
    versions = REPROCESS10X_PARSEMETADATA.out.versions.first()
        .mix(
            REPROCESS10X_LOADDATA.out.versions.first(),
            REPROCESS10X_BAM2FASTQ.out.versions.first(),
            REPROCESS10X_SRA2FASTQ.out.versions.first()
        )

    emit:
    fastq    = fastq_files
    metadata = metadata
    versions = versions
}
