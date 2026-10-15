// IMPORT SUBWORKFLOW
include { REPROCESS10X } from './workflow/main.nf'
include { REPROCESS10X_COLLATE } from './modules/local/reprocess10x/collate'
include { REPROCESS10X_RUNREPORT } from './modules/local/reprocess10x/runreport'

// HELP MESSAGE
def helpMessage() {
    log.info"""
    =============================
    Reprocess Public 10X Datasets
    =============================
    This pipeline loads and reprocesses public 10x datasets from GEO, SRA, ENA, or ArrayExpress.

    Usage: nextflow run main.nf [OPTIONS]

    Required:
        --datasets          Path to a TSV file with dataset and sample IDs.

    Optional:
        --outdir            Directory to save results (default: results).
        --default_specie    Species to assign to samples when metadata is missing or unknown.
                            Accepted values: 'human', 'mouse'. If not set, such samples are skipped
                            by STARsolo with a warning.
        --no_infer_specie   Do not read species from metadata; assign --default_specie to all samples.
                            Requires --default_specie.
        --metaonly          Only fetch metadata, skip downloading and alignment (default: false).
        --starsolo          Run STARsolo alignment after downloading (default: false).
        --cellranger        Run Cell Ranger alignment after downloading (default: false, not yet implemented).
        --help              Show this help message and exit.

    Example:
        nextflow run main.nf \\
          --datasets examples/datasets.tsv \\
          --outdir results \\
          --default_specie human \\
          -resume

    == datasets.tsv format ==
    dataset_id	sample_id
    GSE230685	GSM7232572,GSM7232573
    E-MTAB-9221	ERS4689152,ERS4689153
    PRJEB37166	ERS4605100,ERS4605101

    Leaving sample_id empty processes the WHOLE series. That is intentional, but a
    series is often mixed — use --metaonly first to see what it holds.
    ==========================

    == Output structure ==
    results/
    ├── raw/
    │   └── <dataset_id>/
    │       ├── fastq/
    │       │   └── <sample_id>/    FASTQs (R1, R2, I1)
    │       ├── bam/
    │       │   └── <sample_id>/    10x BAM files
    │       └── sra/
    │           └── <sample_id>/    SRA files
    ├── starsolo/
    │   └── <dataset_id>/           STARsolo count matrices and QC stats
    ├── metadata/
    │   └── <dataset_id>/           Metadata files (links, parsed TSVs, SOFT, etc.)
    ├── index/
    │   ├── fastq.csv               Index of all published FASTQs
    │   ├── bam.csv                 Index of all published BAMs
    │   ├── sra.csv                 Index of all published SRA files
    │   └── starsolo.csv            Index of all STARsolo outputs
    ├── multiqc/                    Run report: run statistics, sample status, chemistry, mapping
    ├── versions.yml                Software versions
    ├── mapping_qc_stats.tsv        Per-sample STARsolo mapping QC statistics
    ├── datasets.tsv                Per-dataset sample, run and species counts
    ├── samples.tsv                 Per-sample resolved species, chemistry and run origins
    ├── run_chemistry.tsv           Per-run chemistry inference and read lengths
    └── skipped_samples.tsv         Samples dropped for want of a species (absent if none)
    ==========================
    """.stripIndent()
}

def unwrapGroupKeys(Map meta) {
    meta.collectEntries { k, v ->
        [k, v instanceof nextflow.extension.GroupKey ? v.getGroupTarget() : v]
    }
}

// Sample metadata as one row of an aligner index: a fixed set of keys in a fixed order.
// The index writes its header from the first record and then each record's own fields in
// order, so a sample missing a key emits a short row under a header that has the column.
// The chemistry keys are exactly that case — they are set upstream only for samples whose
// chemistry run-level inference could settle.
def alignerIndexRow(Map meta) {
    [
        id        : meta.id,
        dataset_id: meta.dataset_id,
        specie    : meta.specie,
        chemistry : meta.chemistry ?: '',  // Cell Ranger --chemistry
        wl        : meta.wl ?: '',         // STARsolo --wl
    ]
}

workflow {
    main:
    // Validate input parameters and show help if needed
    if (params.help || !params.datasets) {
        helpMessage()
        System.exit(params.help ? 0 : 1)
    }

    // Every boolean flag is read as `value.toString().toLowerCase() == 'true'`, so an
    // unrecognised but truthy value (--starsolo 1, yes, or a typo like ture) reads as
    // FALSE. Left unchecked that downloads an entire batch, aligns nothing and exits 0
    // without saying why. Validate all of them up front, before anything is fetched.
    ['metaonly', 'starsolo', 'cellranger', 'no_infer_specie'].each { flag ->
        def value = params[flag]
        if (!(value instanceof Boolean) && !(value?.toString()?.toLowerCase() in ['true', 'false'])) {
            log.error("Invalid value for --${flag}: '${value}'. Expected a boolean (true/false).")
            System.exit(1)
        }
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
    star_human_reference = params.star_human_reference ? channel.value( tuple( [id: "human"], file( params.star_human_reference, checkIfExists: true )) ) : channel.empty()
    star_mouse_reference = params.star_mouse_reference ? channel.value( tuple( [id: "mouse"], file( params.star_mouse_reference, checkIfExists: true )) ) : channel.empty()
    cr_human_reference = params.cr_human_reference ? channel.value( tuple( [id: "human"], file( params.cr_human_reference, checkIfExists: true )) ) : channel.empty()
    cr_mouse_reference = params.cr_mouse_reference ? channel.value( tuple( [id: "mouse"], file( params.cr_mouse_reference, checkIfExists: true )) ) : channel.empty()

    // Define variables
    def metaonlyflag    = params.metaonly && params.metaonly.toString().toLowerCase() == 'true' ? true : false
    def starsoloflag    = params.starsolo && params.starsolo.toString().toLowerCase() == 'true' ? true : false
    def cellrangerflag  = params.cellranger && params.cellranger.toString().toLowerCase() == 'true' ? true : false
    def no_infer_specie = params.no_infer_specie && params.no_infer_specie.toString().toLowerCase() == 'true' ? true : false
    def defaultspecie   = params.default_specie

    // Run main workflow
    REPROCESS10X(
        datasetlist,
        wl_basedir,
        star_human_reference,
        star_mouse_reference,
        cr_human_reference,
        cr_mouse_reference,
        metaonlyflag,
        no_infer_specie,
        defaultspecie,
        starsoloflag,
        cellrangerflag
    )

    // Run-level tables: what the collectFile blocks used to write (versions.yml,
    // mapping_qc_stats.tsv, skipped_samples.tsv) plus datasets/samples/runs/run_chemistry.
    // Every input is collected, so COLLATE runs once, after the last RENAME10XRUN,
    // RENAME10XSAMPLE and MAPPINGQC; ifEmpty([]) keeps it running when a stream is empty
    // (--metaonly, or no aligner requested). Sorted, so the task hash and -resume are stable.
    metadata_files = REPROCESS10X.out.metadata
        .flatMap { _meta, files -> (files instanceof List) ? files : [files] }

    REPROCESS10X_COLLATE(
        datasetlist,
        metadata_files.filter { f -> f.name == 'links.tsv' }.collect(sort: true).ifEmpty([]),
        metadata_files.filter { f -> f.name.endsWith('.sample_x_run.tsv') }.collect(sort: true).ifEmpty([]),
        REPROCESS10X.out.chemistry.map { _meta, json -> json }.collect(sort: true).ifEmpty([]),
        REPROCESS10X.out.soloqc.map { _meta, tsv -> tsv }.collect(sort: true).ifEmpty([]),
        REPROCESS10X.out.versions.collect(sort: true).ifEmpty([]),
        REPROCESS10X.out.fastq.map { meta, _fastqs -> unwrapGroupKeys(meta) }.toSortedList { a, b -> a.id <=> b.id },
        REPROCESS10X.out.skipped.map { meta -> unwrapGroupKeys(meta) }.toSortedList { a, b -> a.id <=> b.id }
    )

    // The run report. It reads this run's own trace, which nextflow.config writes to
    // reports/execution_trace_<trace_report_suffix>.txt (keep the two in step). Its inputs
    // collect COLLATE's tables and every aligner output, so by the time it starts every
    // other task has finished and the trace holds all of them.
    trace_file = file("${workflow.launchDir}/reports/execution_trace_${params.trace_report_suffix}.txt")
    run_info = [
        runName    : workflow.runName,
        sessionId  : workflow.sessionId.toString(),
        commitId   : workflow.commitId ?: '',
        commandLine: workflow.commandLine,
        start      : workflow.start.toString(),
    ]

    REPROCESS10X_RUNREPORT(
        REPROCESS10X_COLLATE.out.tables,
        REPROCESS10X.out.starsolo.map { _meta, dir -> dir }.collect().ifEmpty([]),
        REPROCESS10X.out.cellranger.map { _meta, dir -> dir }.collect().ifEmpty([]),
        trace_file,
        run_info,
        REPROCESS10X.out.planned
            .map { meta -> [meta.id, meta.resources?.reads, meta.resources?.mem_gb?.getAt(0), meta.resources?.hours?.getAt(0)] }
            .collect(flat: false)
            .ifEmpty([])
    )

    publish:
    metadata       = REPROCESS10X.out.metadata.map { meta, files -> tuple(unwrapGroupKeys(meta), files) }
    bam            = REPROCESS10X.out.bam.map { meta, bam -> unwrapGroupKeys(meta) + [path: bam] }
    sra            = REPROCESS10X.out.sra.map { meta, sra -> unwrapGroupKeys(meta) + [path: sra]}
    original_fastq = REPROCESS10X.out.original_fastq.map { meta, fastqs -> unwrapGroupKeys(meta) + [paths: fastqs] }
    runs           = REPROCESS10X.out.runs.map { meta, fastqs -> unwrapGroupKeys(meta) + [paths: fastqs] }
    fastq          = REPROCESS10X.out.fastq
        .flatMap { meta, fastqs -> fastqs.collect { file -> [meta, file] } }
        .map { meta, fastq -> unwrapGroupKeys(meta) + [path: fastq] }
    starsolo   = REPROCESS10X.out.starsolo.map { meta, starsolo -> alignerIndexRow(unwrapGroupKeys(meta)) + [path: starsolo] }
    soloqc     = REPROCESS10X.out.soloqc.map { meta, soloqc -> meta.getGroupTarget() + [path: soloqc] }
    cellranger = REPROCESS10X.out.cellranger.map { meta, cellranger -> alignerIndexRow(unwrapGroupKeys(meta)) + [path: cellranger] }
    tables     = REPROCESS10X_COLLATE.out.tables.flatten()
    report     = REPROCESS10X_RUNREPORT.out.report.mix(REPROCESS10X_RUNREPORT.out.data)
}

output {
    metadata {
        label "metadata"
        path { meta, _files -> "metadata/${meta.id}" }
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
            path "index/renamed_fastq.csv"
            header true
            sep ','
        }
        path { output -> "raw/${output.dataset_id}/renamed_fastq/${output.id}" }
    }
    original_fastq {
        label "fastq"
        label "raw"
        index {
            path "index/original_fastq.csv"
            header true
            sep ','
        }
        path { output -> "raw/${output.dataset_id}/fastq/${output.sample_id}/${output.id}" }
    }
    runs {
        label "fastq"
        label "raw"
        index {
            path "index/runs.csv"
            header true
            sep ','
        }
        path { output -> "raw/${output.dataset_id}/runs/${output.sample_id}/${output.id}" }
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
    cellranger {
        label "cellranger"
        index {
            path "index/cellranger.csv"
            header true
            sep ','
        }
        path { output -> "cellranger/${output.dataset_id}/" }
    }
    tables {
        // at the outdir root, where mapping_qc_stats.tsv & co. have always been; a plain
        // `path "."` would keep COLLATE's tables/ subdir
        path { table -> table >> table.name }
    }
    report {
        path "multiqc"
    }

}