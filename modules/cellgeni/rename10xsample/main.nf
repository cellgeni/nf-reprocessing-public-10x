/*
 * Module: cellgeni/rename10xsample
 */

process RENAME10XSAMPLE {
    tag "${meta.id}"
    container "quay.io/cellgeni/metacells-python:latest"

    input:
    // run_jsons are the per-run chemistry reports from RENAME10XRUN. They let this
    // step check CB/UMI geometry across a sample's runs instead of re-deriving read
    // lengths from a small head sample, which is what made runs of one sample look
    // like they disagreed. Runs that never went through RENAME10XRUN (BAM-derived
    // ones) simply contribute no JSON, so the list may be empty or partial.
    tuple val(meta), path(fastqs, stageAs: "fastqs/*"), path(run_jsons, stageAs: "run_jsons/*")

    output:
    tuple val(meta), path("*_R*_001.fastq.gz"), emit: reads
    tuple val(meta), path("*_I*_001.fastq.gz"), optional: true, emit: index
    path  "versions.yml",                               emit: versions

    script:
    def args = task.ext.args ?: ""
    def run_jsons_arg = run_jsons ? "--run-jsons ${run_jsons}" : ""
    """
    rename_fastqs_recommended.py \\
        ${args} \\
        --fastqs ${fastqs} \\
        ${run_jsons_arg} \\
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
