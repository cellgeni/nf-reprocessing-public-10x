/*
 * Module: cellgeni/rename10xsample
 */

process RENAME10XSAMPLE {
    tag "${meta.id}"
    container "quay.io/cellgeni/metacells-python:latest"

    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")

    output:
    tuple val(meta), path("*_R*_001.fastq.gz"), emit: reads
    tuple val(meta), path("*_I*_001.fastq.gz"), optional: true, emit: index
    path  "versions.yml",                               emit: versions

    script:
    def args = task.ext.args ?: ""
    """
    rename_fastqs.py \\
        ${args} \\
        --fastqs ${fastqs} \\
        --sample-id ${meta.id} \\
        --outdir .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo "stub" | gzip > ${meta.id}_S1_L1_R1_001.fastq.gz
    echo "stub" | gzip > ${meta.id}_S1_L1_R2_001.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
