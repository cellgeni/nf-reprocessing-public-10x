include { REPROCESS10X_PARSEMETADATA } from '../../../modules/local/reprocess10x/parsemetadata'
include { REPROCESS10X_LOADDATA } from '../../../modules/local/reprocess10x/loaddata'

workflow DOWNLOAD10X {

    take:
    datasets     // channel: [ val(meta), [ sample_ids ] ]
    wl_basedir   // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)
    
    main:
    
    REPROCESS10X_PARSEMETADATA(datasets)

    // Get links from metadata
    links = REPROCESS10X_PARSEMETADATA.out.links
                                          .splitCsv(sep: '\t', strip: true)
                                          .map { meta, row -> 
                                              [
                                                [
                                                  id       : row[0],
                                                  sample_id: row[4],
                                                  dataset_id : meta.id,
                                                  specie   : row[1],
                                                  type     : row[3]
                                                ],
                                                row[2].split(";")
                                            ]
                                          }

      
    links.view()
    // Load data from links
    REPROCESS10X_LOADDATA(links)

    REPROCESS10X_LOADDATA.out.urls.view()
    REPROCESS10X_LOADDATA.out.fastq.view()
    REPROCESS10X_LOADDATA.out.sra.view()

    // // Convert data if needed
    // REPROCESS10X_BAM2FASTQ(REPROCESS10X_LOADDATA.out.bam)
    // REPROCESS10X_SRA2FASTQ(REPROCESS10X_BAM2FASTQ.out.sra)

    // // Combine all fastq files
    // fastq_files = REPROCESS10X_LOADDATA.out.fastq
    //                                        .concat(
    //                                               REPROCESS10X_BAM2FASTQ.out.fastq,
    //                                               REPROCESS10X_SRA2FASTQ.out.fastq
    //                                        )
    // // Rename fastq files according to Cell Ranger format
    // REPROCESS10X_RENAMEFASTQ(fastq_files)

    // // Group fastq files by sample
    // sample_fastqs = REPROCESS10X_RENAMEFASTQ.fastq.map { run_meta, fastq -> 
    //                                     [
    //                                         id: run_meta.sample_id,
    //                                         specie: run_meta.specie,
    //                                     ]
    // }


    // emit:
    // fastq    = sample_fastqs                      // channel: [ val(meta), [ bam ] ]
    // tsv      = REPROCESS10X_LOADMETADATA.out.tsv  // channel: [ val(meta), [ tsv ] ]
    // versions = ch_versions                        // channel: [ versions.yml ]
}
