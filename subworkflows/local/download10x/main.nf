include { REPROCESS10X_PARSEMETADATA } from '../../../modules/local/reprocess10x/parsemetadata'
include { REPROCESS10X_LOADDATA } from '../../../modules/local/reprocess10x/loaddata'
include { REPROCESS10X_BAM2FASTQ } from '../../../modules/local/reprocess10x/bam2fastq'
include { REPROCESS10X_SRA2FASTQ } from '../../../modules/local/reprocess10x/sra2fastq'

workflow DOWNLOAD10X {

    take:
    datasets     // channel: [ val(meta), [ sample_ids ] ]
    wl_basedir   // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)
    
    main:
    
    REPROCESS10X_PARSEMETADATA(datasets)

    // Get links for each run file from metadata
    links = REPROCESS10X_PARSEMETADATA.out.links
                                          .splitCsv(sep: '\t', strip: true)
                                          .map { meta, row -> 
                                              [
                                                [
                                                  id       : row[0], // run ID
                                                  sample_id: row[4], // sample ID
                                                  dataset_id : meta.id, // dataset's series ID
                                                  specie   : row[1], // specie
                                                  type     : row[3] // file type i.e. BAM, FASTQ, SRA ...
                                                ],
                                                row[2].split(";") // split URLs by semicolon in case of multiple .fastq files
                                            ]
                                          }

      
    links.view()
    // Load data from links
    REPROCESS10X_LOADDATA(links)

    REPROCESS10X_LOADDATA.out.urls.view()
    REPROCESS10X_LOADDATA.out.fastq.view()
    REPROCESS10X_LOADDATA.out.sra.view()
    REPROCESS10X_LOADDATA.out.bam.view()

    // // Convert data if needed
    REPROCESS10X_BAM2FASTQ(REPROCESS10X_LOADDATA.out.bam)
    REPROCESS10X_SRA2FASTQ(REPROCESS10X_LOADDATA.out.sra, wl_basedir)

    // // Combine all fastq channels and group by sample
    fastq_files = REPROCESS10X_LOADDATA.out.fastq
                          .mix(
                                REPROCESS10X_BAM2FASTQ.out.fastq,
                                REPROCESS10X_SRA2FASTQ.out.fastq
                          )
                          // Leave only sample id and dataset id in metadata
                          .map { run_meta, fastq ->
                              [
                                [
                                  id       : run_meta.sample_id,
                                  dataset_id : run_meta.dataset_id,
                                ],
                                fastq
                              ]
                          }
                          // Group by sample id and dataset id
                          .groupTuple()
    fastq_files.view()
    // // Rename fastq files according to Cell Ranger format
    // REPROCESS10X_RENAMEFASTQ(fastq_files)


    // emit:
    // fastq    = REPROCESS10X_RENAMEFASTQ.fastq     // channel: [ val(meta), [ bam ] ]
    // tsv      = REPROCESS10X_LOADMETADATA.out.tsv  // channel: [ val(meta), [ tsv ] ]
    // versions = ch_versions                        // channel: [ versions.yml ]
}
