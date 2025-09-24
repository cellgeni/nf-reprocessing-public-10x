// IMPORT SUBWORKFLOW
include { DOWNLOAD10X } from './subworkflows/local/download10x/'
include { STARSOLO10X as STARSOLO10X_HUMAN } from './subworkflows/local/starsolo10x/'
include { STARSOLO10X as STARSOLO10X_MOUSE } from './subworkflows/local/starsolo10x/'


// HELP MESSAGE
def helpMessage() {
    log.info"""
    =============================
    Reprocess Public 10X Datasets
    =============================
    This pipeline performs loading and processing of public datasets stored on GEO, SRA, ENA or ArrayExpress.
    Usage: nextflow run main.nf [OPTIONS]
        options:
            --datasets   Path to a .tsv file containing dataset information (e.g., datasets.tsv).
                         First column should be 'dataset_id' and second column should be a list 
                         of comma separated sample IDs.
            --mapper     Specify the mapper to use (e.g., 'starsolo' or 'cellranger').
            --help       Show this help message and exit.

    Examples:
        nextflow run main.nf --datasets examples/datasets.tsv
    == samples.csv format ==
    dataset_id	sample_id
    GSE230685	GSM7232572,GSM7232573
    E-MTAB-9221	ERS4689152,ERS4689153
    ========================
    """.stripIndent()
}

workflow {
    if (params.help) {
        helpMessage()
        System.exit(0)
    } 
    // STEP 0: Validate input parameters
    // Convert dataset list to channels
    datasets = Channel
        .fromPath(params.datasets)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_list =  row.sample_id.split(',')
            [
                [id: groupKey(row.dataset_id, sample_list.size())],
                sample_list
            ]
        }
    
    // STEP 1: Download datasets
    DOWNLOAD10X(
        datasets,
        params.wl_basedir
    )

    // Group samples by specie
    fastqs = DOWNLOAD10X.out.fastq
        .branch { meta, _fastqs ->
            human: meta.specie == 'Homo sapiens'
            mouse: meta.specie == 'Mus musculus'
            other: true
        }

    //fastqs.human.view { meta, fastq -> "HUMAN: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastq (${fastq.getClass().simpleName})" }
    //fastqs.mouse.view { meta, fastq -> "MOUSE: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastq (${fastq.getClass().simpleName})" }
    //fastqs.other.view { meta, fastq -> "OTHER: meta=[${meta.collect { k, v -> "$k: $v (${v.getClass().simpleName})" }.join(', ')}], fastqs=$fastq (${fastq.getClass().simpleName})" }

    // Warn user if there are unexpected species detected
    fastqs.other.count().subscribe { count ->
        if (count > 0) {
            log.warn "Detected ${count} fastq files from unexpected species"
        }
    }
    
    // STEP2: Run STARsolo on fastq files for human and mouse samples
    human_reference = channel.value(tuple( [id: "human"], file( params.human_reference ) ))
    mouse_reference = channel.value(tuple( [id: "mouse"], file( params.mouse_reference ) ))

    STARSOLO10X_HUMAN(fastqs.human, human_reference, params.wl_basedir)
    STARSOLO10X_MOUSE(fastqs.mouse, mouse_reference, params.wl_basedir)

    // STEP 3: Collect outputs
    // Collect mapping QC
    STARSOLO10X_HUMAN.out.qc_stats
        .mix(STARSOLO10X_MOUSE.out.qc_stats)
        .splitCsv(sep: '\t', skip: 1)
        .collectFile(
            name: 'mapping_qc_stats.tsv',
            storeDir: params.output_dir,
            newLine: true,
            sort: { line -> line.split('\t').take(2).join('') },  // sort by first two columns (Dataset, sample ID),
            seed: "Dataset\tSample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"
        ) { meta, row -> 
            "${meta.id}\t${row.join('\t')}"
        }
        .subscribe { __ -> 
                log.info("Mapping QC stats saved to ${params.output_dir}/mapping_qc_stats.tsv")
            }

    // Collect versions
    DOWNLOAD10X.out.versions
        .mix(
            STARSOLO10X_HUMAN.out.versions,
            STARSOLO10X_MOUSE.out.versions
        )
        .splitText(by: 20)
        .unique()
        .collectFile(name: 'versions.yml', storeDir: params.output_dir, sort: true)
        .subscribe { __ -> 
                log.info("Versions saved to ${params.output_dir}/versions.yml")
            }
}