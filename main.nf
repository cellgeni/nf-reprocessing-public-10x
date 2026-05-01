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
    if (params.help || !params.datasets) {
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

    // Resolve effective specie for each sample:
    //   --specie human/mouse  → override all samples to that specie
    //   --specie auto         → use meta.specie from metadata; fall back to --default_specie if blank/NULL/UNKNOWN
    def unknown_values = [null, '', 'NULL', 'UNKNOWN']
    resolved_fastqs = DOWNLOAD10X.out.fastq
        .map { meta, fastqs ->
            def effective_specie
            if (params.specie == 'human') {
                effective_specie = 'Homo sapiens'
            } else if (params.specie == 'mouse') {
                effective_specie = 'Mus musculus'
            } else {
                effective_specie = (meta.specie in unknown_values) ? params.default_specie : meta.specie
            }
            if (effective_specie in unknown_values) {
                log.warn "Sample ${meta.id} (${meta.dataset_id}) has unknown species and no --default_specie set — skipping STARsolo"
            }
            [meta + [specie: effective_specie], fastqs]
        }

    // Group samples by specie
    fastqs = resolved_fastqs
        .branch { meta, _fastqs ->
            human: meta.specie == 'Homo sapiens'
            mouse: meta.specie == 'Mus musculus'
            other: true
        }

    fastqs.other
        .map { meta, _fastqs ->
            log.warn "Sample ${meta.id} (${meta.dataset_id}) has unexpected species '${meta.specie}' — skipping STARsolo"
            [meta, _fastqs]
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