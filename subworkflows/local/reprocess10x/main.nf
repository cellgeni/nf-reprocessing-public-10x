include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow REPROCESS10X {

    take:
    datasets // channel: [ val(meta), [ sample_ids ] ]

    main:
    
    REPROCESS10X_LOADMETADATA(datasets)

    // Get links from metadata
    links = REPROCESS10X_LOADMETADATA.out.links
                                          .splitCsv(sep: '\t', strip: true)
                                          .map { row -> 
                                              [
                                                [
                                                  id: row[1],
                                                  sample_id: row[0],
                                                  specie: row[2],
                                                  type: row[4]
                                                ],
                                                row[3]
                                            ]
                                          }
    
    // Load data from links
    REPROCESS10X_LOADDATA(links)

    // Convert data if needed
    REPROCESS10X_BAM2FASTQ(REPROCESS10X_LOADDATA.out.bam)
    REPROCESS10X_SRA2FASTQ(REPROCESS10X_BAM2FASTQ.out.sra)

    // Combine all fastq files
    fastq_files = REPROCESS10X_LOADDATA.out.fastq
                                           .concat(
                                                  REPROCESS10X_BAM2FASTQ.out.fastq,
                                                  REPROCESS10X_SRA2FASTQ.out.fastq
                                           )
    
    // Group fastq files by sample
    sample_fastqs = fastq_files.map { meta, fastq ->
        [
            meta.sample_id,
            meta,
            fastq
        ]
    }

    
                        



    

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
