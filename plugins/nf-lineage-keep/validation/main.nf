/*
 * End-to-end check of nf-lineage-keep against the pipeline's shape, in miniature and
 * with no network: per-dataset metadata -> per-run download -> per-run rename ->
 * per-sample merge -> per-sample align -> per-dataset QC, plus two end-of-run barriers,
 * publishing every stage by hard link as main.nf does. ALIGN fails for good on
 * --fail (after one retry, as the pipeline's errorStrategy does) and leaves scratch
 * behind on every attempt, like STAR's _STARtmp, plus a partial BAM. DOWNLOAD of the
 * runs in --flaky fails once and succeeds on its retry.
 *
 *   ./run.sh            # see run.sh: a run, checks, then a -resume of rerun.tsv
 */

// two samples per dataset, two runs per sample
def layout() {
    return [GSE1: ['S1', 'S2'], GSE2: ['S3', 'S4']]
}

process FETCH {
    tag "${meta.id}"
    input:
    tuple val(meta), val(sample_ids)
    output:
    tuple val(meta), path('links.tsv'), emit: links
    path 'versions.yml', emit: versions
    script:
    def wanted = sample_ids ? sample_ids.tokenize(',') : layout()[meta.id.toString()]
    def rows = wanted.collectMany { s -> ['a', 'b'].collect { r -> "${s}${r}\t${s}" } }.join('\n')
    """
    printf '${rows}\\n' > links.tsv
    echo 'FETCH: 1' > versions.yml
    """
}

// a run in --flaky fails its first attempt, leaving a partial download and a small
// report behind, and succeeds on the retry
process DOWNLOAD {
    tag "${meta.id}"
    errorStrategy { task.attempt == 1 ? 'retry' : 'ignore' }
    input:
    val meta
    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: fastq
    path 'versions.yml', emit: versions
    script:
    def flaky = task.attempt == 1 && params.flaky.tokenize(',').contains(meta.id)
    """
    if ${flaky}; then
        head -c 1500000 /dev/urandom > partial.fastq.gz
        echo 'connection reset' > report.txt
        echo "simulated download failure" >&2; exit 1
    fi
    head -c 200000 /dev/urandom > ${meta.id}.fastq.gz
    echo 'DOWNLOAD: 1' > versions.yml
    """
}

process RENAME_RUN {
    tag "${meta.id}"
    input:
    tuple val(meta), path(fq, stageAs: 'in/*')
    output:
    tuple val(meta), path("${meta.id}_R1_001.fastq.gz"), emit: reads
    path 'versions.yml', emit: versions
    script:
    // a hard link, as infer_10x_run.py makes by default
    """
    ln \$(readlink -f ${fq}) ${meta.id}_R1_001.fastq.gz
    echo 'RENAME_RUN: 1' > versions.yml
    """
}

process RENAME_SAMPLE {
    tag "${meta.id}"
    input:
    tuple val(meta), path(fqs, stageAs: 'in/*')
    output:
    tuple val(meta), path("${meta.id}_S1_R1_001.fastq.gz"), emit: reads
    path 'versions.yml', emit: versions
    script:
    """
    cat ${fqs} > ${meta.id}_S1_R1_001.fastq.gz
    echo 'RENAME_SAMPLE: 1' > versions.yml
    """
}

process ALIGN {
    tag "${meta.id}"
    errorStrategy { task.attempt == 1 ? 'retry' : 'ignore' }
    input:
    tuple val(meta), path(fq, stageAs: 'fastqs/*')
    output:
    tuple val(meta), path("${meta.id}"), emit: mapping
    path 'versions.yml', emit: versions
    script:
    """
    mkdir -p ${meta.id}/_STARtmp
    head -c 300000 /dev/urandom > ${meta.id}/_STARtmp/scratch.bin
    if [[ "${meta.id}" == "${params.fail}" ]]; then
        head -c 1500000 /dev/urandom > ${meta.id}/partial.bam
        echo 'EXITING because of fatal error' > ${meta.id}/Log.out
        echo "simulated failure" >&2; exit 1
    fi
    wc -c < fastqs/* > ${meta.id}/matrix.mtx
    rm -rf ${meta.id}/_STARtmp
    echo 'ALIGN: 1' > versions.yml
    """
}

process QC {
    tag "${meta.id}"
    input:
    tuple val(meta), path(samples)
    output:
    tuple val(meta), path("${meta.id}.qc.tsv"), emit: tsv
    script:
    """
    ls ${samples} > ${meta.id}.qc.tsv
    """
}

process COLLATE {
    input:
    path versions, stageAs: 'v*/versions.yml'
    path qc
    output:
    path 'versions_all.yml', emit: tables
    script:
    """
    cat ${versions} > versions_all.yml
    """
}

process REPORT {
    input:
    path mappings
    output:
    path 'report.txt', emit: report
    script:
    """
    ls -d ${mappings} > report.txt
    """
}

workflow {
    main:
    datasets = channel.fromPath(params.datasets)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [[id: row.dataset_id], row.sample_id ?: ''] }

    FETCH(datasets)

    runs = FETCH.out.links.flatMap { meta, links ->
        def rows = links.splitCsv(sep: '\t')
        def bySample = rows.groupBy { r -> r[1] }
        rows.collect { r ->
            [id: r[0], sample_id: groupKey(r[1], bySample[r[1]].size()), dataset_id: groupKey(meta.id, bySample.size())]
        }
    }
    DOWNLOAD(runs)
    RENAME_RUN(DOWNLOAD.out.fastq)

    samples = RENAME_RUN.out.reads
        .map { m, fq -> [groupKey([id: m.sample_id.getGroupTarget(), dataset_id: m.dataset_id], m.sample_id.getGroupSize()), fq] }
        .groupTuple(remainder: true)
        .map { k, fqs -> [k.getGroupTarget(), fqs] }
    RENAME_SAMPLE(samples)

    ALIGN(RENAME_SAMPLE.out.reads)

    byDataset = ALIGN.out.mapping
        .map { m, dir -> [groupKey([id: m.dataset_id.getGroupTarget()], m.dataset_id.getGroupSize()), dir] }
        .groupTuple(remainder: true)
    QC(byDataset)

    versions = FETCH.out.versions.first()
        .mix(DOWNLOAD.out.versions.first(), RENAME_RUN.out.versions.first(),
             RENAME_SAMPLE.out.versions.first(), ALIGN.out.versions.first())
    COLLATE(versions.collect(), QC.out.tsv.map { _m, f -> f }.collect().ifEmpty([]))
    REPORT(ALIGN.out.mapping.map { _m, d -> d }.collect().ifEmpty([]))

    publish:
    raw     = DOWNLOAD.out.fastq.map { m, f -> [dataset: m.dataset_id.toString(), sample: m.sample_id.toString(), path: f] }
    runs    = RENAME_RUN.out.reads.map { m, f -> [dataset: m.dataset_id.toString(), sample: m.sample_id.toString(), path: f] }
    fastq   = RENAME_SAMPLE.out.reads.map { m, f -> [dataset: m.dataset_id.toString(), sample: m.id, path: f] }
    mapping = ALIGN.out.mapping.map { m, d -> [dataset: m.dataset_id.toString(), path: d] }
    qc      = QC.out.tsv.map { m, f -> [dataset: m.id.toString(), path: f] }
    tables  = COLLATE.out.tables
    report  = REPORT.out.report
}

output {
    raw     { path { o -> "raw/${o.dataset}/fastq/${o.sample}" } }
    runs    { path { o -> "raw/${o.dataset}/runs/${o.sample}" } }
    fastq   { path { o -> "raw/${o.dataset}/renamed_fastq/${o.sample}" } }
    mapping { path { o -> "starsolo/${o.dataset}" } }
    qc      { path { o -> "starsolo/${o.dataset}" } }
    tables  { path '.' }
    report  { path '.' }
}
