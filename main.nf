// IMPORT SUBWORKFLOW
include { DOWNLOAD10X } from './subworkflows/local/download10x/'
include { STARSOLO10X as STARSOLO10X_HUMAN } from './subworkflows/local/starsolo10x/'
include { STARSOLO10X as STARSOLO10X_MOUSE } from './subworkflows/local/starsolo10x/'
include { REPROCESS10X } from './workflow/main.nf'

// HELP MESSAGE
def helpMessage() {
    log.info"""
    =============================
    Reprocess Public 10X Datasets
    =============================
    This pipeline performs loading and processing of public datasets stored on GEO, SRA, ENA or ArrayExpress.

    Usage: nextflow run main.nf [OPTIONS]

    Required:
        --datasets          Path to a TSV file with dataset and sample IDs.

    Optional:
        --specie            Species to use for all samples: 'human', 'mouse', or 'auto' (default: auto).
                            'auto' reads species from metadata; falls back to --default_specie if unknown.
        --default_specie    Fallback species when metadata is missing or unknown (e.g. "human" or "mouse").
        --output_dir        Directory to save results (default: results).
        --help              Show this help message and exit.

    Example:
        nextflow run main.nf \\
          -profile standard \\
          --datasets examples/datasets.tsv \\
          --output_dir results \\
          -resume

    == datasets.tsv format ==
    dataset_id	sample_id
    GSE230685	GSM7232572,GSM7232573
    E-MTAB-9221	ERS4689152,ERS4689153
    PRJEB37166	ERS4605100,ERS4605101
    ==========================
    """.stripIndent()
}

def unwrapGroupKeys(Map meta) {
    meta.collectEntries { k, v ->
        [k, v instanceof nextflow.extension.GroupKey ? v.getGroupTarget() : v]
    }
}

workflow {
    main:
    // Validate input parameters and show help if needed
    if (params.help || !params.datasets) {
        helpMessage()
        System.exit(params.help ? 0 : 1)
    }

    if (!(params.metaonly instanceof Boolean)) {
        log.error("Invalid value for --metaonly: ${params.metaonly}. Expected a boolean (true/false).")
        System.exit(1)
    }

    if (params.no_infer_specie && params.default_specie == null) {
        log.error("When --no_infer_specie is set, --default_specie must be provided.")
        System.exit(1)
    }

    if ( params.default_specie && !['human', 'mouse'].contains(params.default_specie.toLowerCase()) ) {
        log.error("Invalid value for --default_specie: ${params.default_specie}. Expected 'human' or 'mouse'.")
        System.exit(1)
    }

    // Load files
    datasetlist     = channel.value( file( params.datasets, checkIfExists: true ) )
    wl_basedir      = params.wl_basedir ? channel.value( file( params.wl_basedir, checkIfExists: true ) ) : channel.empty()
    human_reference = params.human_reference ? channel.value( tuple( [id: "human"], file( params.human_reference, checkIfExists: true )) ) : channel.empty()
    mouse_reference = params.mouse_reference ? channel.value( tuple( [id: "mouse"], file( params.mouse_reference, checkIfExists: true )) ) : channel.empty()

    // Define variables
    def metaonlyflag    = params.metaonly ? true : false
    def no_infer_specie = params.no_infer_specie ? false : true
    def defaultspecie   = params.default_specie

    // Run main workflow
    REPROCESS10X(
        datasetlist,
        wl_basedir,
        human_reference,
        mouse_reference,
        metaonlyflag,
        no_infer_specie,
        defaultspecie
    )

    // Collect versions
    REPROCESS10X.out.versions
        .splitText(by: 20)
        .unique()
        .collectFile(name: 'versions.yml', storeDir: params.outdir, sort: true)
        .subscribe { __ -> 
                log.info("Versions saved to ${params.outdir}/versions.yml")
            }
    
    // Collect mapping QC stats
    REPROCESS10X.out.soloqc
        .splitCsv(sep: '\t', skip: 1)
        .collectFile(
            name: 'mapping_qc_stats.tsv',
            storeDir: params.outdir,
            newLine: true,
            seed: "Dataset\tSample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"
        ) { meta, row -> 
            "${meta.id}\t${row.join('\t')}"
        }
        .subscribe { __ -> 
                log.info("Mapping QC stats saved to ${params.outdir}/mapping_qc_stats.tsv")
            }
    
    publish:
    metadata = REPROCESS10X.out.metadata.map { meta, files -> tuple(unwrapGroupKeys(meta), files) }
    bam      = REPROCESS10X.out.bam.map { meta, bam -> unwrapGroupKeys(meta) + [path: bam] }
    sra      = REPROCESS10X.out.sra.map { meta, sra -> unwrapGroupKeys(meta) + [path: sra]}
    fastq    = REPROCESS10X.out.fastq
        .flatMap { meta, fastqs -> fastqs.collect { file -> [meta, file] } }
        .map { meta, fastq -> unwrapGroupKeys(meta) + [path: fastq] }
    starsolo = REPROCESS10X.out.starsolo.map { meta, starsolo -> unwrapGroupKeys(meta) + [path: starsolo] }
    soloqc   = REPROCESS10X.out.soloqc.map { meta, soloqc -> meta.getGroupTarget() + [path: soloqc] }
}

output {
    metadata {
        label "metadata"
        path { meta, files -> "metadata/${meta.id}" }
    }
    bam {
        label "bam"
        label "raw"
        index {
            path "index/bam.csv"
            header true
            sep ','
        }
        path { output -> "raw/${output.dataset_id}/bam/${output.sample_id}" }
    }
    sra {
        label "sra"
        label "raw"
        index {
            path "index/sra.csv"
            header true
            sep ','
        }
        path { output -> "raw/${output.dataset_id}/sra/${output.sample_id}" }
    }
    fastq {
        label "fastq"
        label "raw"
        index {
            path "index/fastq.csv"
            header true
            sep ','
        }
        path { output -> "raw/${output.dataset_id}/fastq/${output.id}" }
    }
    starsolo {
        label "starsolo"
        index {
            path "index/starsolo.csv"
            header true
            sep ','
        }
        path { output -> "starsolo/${output.dataset_id}/" }
    }
    soloqc {
        label "metadata"
        label "qc"
        index {
            path "index/qc.csv"
            header true
            sep ','
        }
        path { output -> "starsolo/${output.id}/" }
    }

}