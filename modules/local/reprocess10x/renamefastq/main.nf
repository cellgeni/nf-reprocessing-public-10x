process REPROCESS10X_RENAMEFASTQ {
    tag "Inferring 10x FASTQ layout and renaming files for $meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x:latest':
        'quay.io/cellgeni/reprocess_10x:latest' }"

    input:
    tuple val(meta), path(fastqs)
    path whitelists

    output:
    tuple val(meta), path("fastqs/*.fastq.gz"), emit: fastq
    tuple val(meta), path("10x_fastq_inference.json"), emit: report
    tuple val(meta), path("10x_fastq_inference.tsv"),  emit: summary
    path "versions.yml",                            emit: versions

    script:
    def fq_args = fastqs instanceof List ? fastqs.collect { "\"${it}\"" }.join(' ') : "\"${fastqs}\""
    """
    mkdir -p fastqs

    infer_10x_fastqs.py \
        --sample "${meta.id}" \
        --whitelist-dir "${whitelists}" \
        --outdir fastqs \
        --target "${params.fastq_rename_target}" \
        --json-report 10x_fastq_inference.json \
        --tsv-report 10x_fastq_inference.tsv \
        ${fq_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{ print \$2 }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p fastqs
    touch fastqs/${meta.id}_S1_L001_R1_001.fastq.gz
    touch fastqs/${meta.id}_S1_L001_R2_001.fastq.gz
    cat <<-END_JSON > 10x_fastq_inference.json
    {"sample": "${meta.id}", "stub": true}
    END_JSON
    cat <<-END_TSV > 10x_fastq_inference.tsv
    sample	selected_chemistry	ambiguous_chemistry	safe_rename	warnings
    ${meta.id}	STUB	false	true	
    END_TSV
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
