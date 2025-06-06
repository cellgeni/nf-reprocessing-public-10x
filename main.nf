// IMPORT SUBWORKFLOW
include { DOWNLOAD10X } from './subworkflows/local/download10x/'

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

    // Convert dataset list to channels
    datasets = Channel
        .fromPath(params.datasets)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            [
                [id: row.dataset_id],
                row.sample_id.split(',')
            ]
        }
    
    DOWNLOAD10X(
        datasets,
        params.wl_basedir
    )
}