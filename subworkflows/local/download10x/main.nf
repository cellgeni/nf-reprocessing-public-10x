include { REPROCESS10X_PARSEMETADATA } from '../../../modules/local/reprocess10x/parsemetadata'
include { REPROCESS10X_LOADDATA } from '../../../modules/local/reprocess10x/loaddata'
include { REPROCESS10X_BAM2FASTQ } from '../../../modules/local/reprocess10x/bam2fastq'
include { REPROCESS10X_SRA2FASTQ } from '../../../modules/local/reprocess10x/sra2fastq'
include { REPROCESS10X_STARSOLO as REPROCESS10X_STARSOLO_HUMAN } from '../../../modules/local/reprocess10x/starsolo'
include { REPROCESS10X_STARSOLO as REPROCESS10X_STARSOLO_MOUSE } from '../../../modules/local/reprocess10x/starsolo'

workflow DOWNLOAD10X {

    take:
    datasets     // channel: [ val(meta), [ sample_ids ] ]
    wl_basedir   // channel: [ dirpath ] a path to whitelist base directory (see https://github.com/cellgeni/nf-reprocessing-public-10x/tree/main/data/whitelists/)
    
    main:
    
    REPROCESS10X_PARSEMETADATA(datasets)

    // Get links for each run file from metadata
    links = REPROCESS10X_PARSEMETADATA.out.links
                                          .splitCsv(sep: '\t', strip: true)
                                          .flatMap { meta, row -> 
                                              def run_meta = [
                                                  id       : row[0], // run ID
                                                  sample_id: row[4], // sample ID
                                                  dataset_id : meta.id, // dataset's series ID
                                                  specie   : row[1], // specie
                                                  type     : row[3] // file type i.e. BAM, FASTQ, SRA ...
                                              ]
                                              def urls = row[2].split(";")
                                              return urls.collect { url -> [run_meta, url] }
                                          }
                                          

      
    links.view()
    // Load data from links
    REPROCESS10X_LOADDATA(links)

    REPROCESS10X_LOADDATA.out.fastq.view()
    REPROCESS10X_LOADDATA.out.sra.view()
    REPROCESS10X_LOADDATA.out.bam.view()

    // Convert data if needed
    REPROCESS10X_BAM2FASTQ(REPROCESS10X_LOADDATA.out.bam)
    REPROCESS10X_SRA2FASTQ(REPROCESS10X_LOADDATA.out.sra, wl_basedir)

    // Combine all fastq channels and group by sample
    fastq_files = REPROCESS10X_LOADDATA.out.fastq
                          // Combine fastq files for each run as they were loaded separately
                          .groupTuple()
                          // Combine fastq files from BAM and SRA conversion
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
                                  specie: run_meta.specie,
                                  type: run_meta.type
                                ],
                                fastq
                              ]
                          }
                          // Group by sample id and dataset id
                          .groupTuple()
                          // Split into separate channels for each specie
                          .branch { meta, fastqs ->
                              human: meta.specie == 'Homo sapiens'
                              mouse: meta.specie == 'Mus musculus'
                              other: true
                          }
    fastq_files.human.view()
    fastq_files.mouse.view()
    fastq_files.other.view()

    // Run STARsolo on fastq files for human and mouse samples
    human_reference = channel.value(file( params.human_reference ))
    mouse_reference = channel.value(file( params.mouse_reference ))

    REPROCESS10X_STARSOLO_HUMAN(fastq_files.human, human_reference)
    REPROCESS10X_STARSOLO_MOUSE(fastq_files.mouse, mouse_reference)

    // Warn user if there are unexpected species detected
    fastq_files.other.count().subscribe { count ->
        if (count > 0) {
            log.warn "Detected ${count} fastq files from unexpected species"
        }
    }

    // emit:
    // fastq    = REPROCESS10X_RENAMEFASTQ.fastq     // channel: [ val(meta), [ bam ] ]
    // tsv      = REPROCESS10X_LOADMETADATA.out.tsv  // channel: [ val(meta), [ tsv ] ]
    // versions = ch_versions                        // channel: [ versions.yml ]
}
