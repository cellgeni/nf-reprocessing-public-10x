/*
 * Module: cellgeni/rename10xrun
 */

process RENAME10XRUN {
    tag "${meta.id}"
    container "quay.io/cellgeni/metacells-python:latest"

    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    path  whitelist_dir

    output:
    tuple val(meta), path("*_R*_001.fastq.gz"), emit: reads
    tuple val(meta), path("*_I*_001.fastq.gz"), optional: true, emit: index
    tuple val(meta), path("${meta.id}.chemistry.json"), emit: chemistry
    path  "versions.yml",                               emit: versions

    script:
    def args = task.ext.args ?: ""
    """
    infer_10x_run_recommended.py \\
        --fastqs ${fastqs} \\
        --run-id ${meta.id} \\
        --whitelist-dir ${whitelist_dir} \\
        ${args} \\
        --outdir . \\
        --json ${meta.id}.chemistry.json \\
        --tsv ${meta.id}.chemistry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def lane = meta.lane ?: 1
    """
    echo "stub" | gzip > ${meta.id}_S1_L1_R1_001.fastq.gz
    echo "stub" | gzip > ${meta.id}_S1_L1_R2_001.fastq.gz
    echo '{"sample":"${meta.id}","chemistry":"stub"}' > ${meta.id}.chemistry.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
