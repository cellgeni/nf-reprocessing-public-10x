process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'
    container "quay.io/nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    tuple val(refmeta), path(reference)

    output:
    tuple val(meta), path("${meta.id}"), emit: mapping
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cellranger count \
        --id="${meta.id}" \
        --fastqs="fastqs" \
        --transcriptome="${reference}" \
        --sample="${meta.id}" \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}/outs/"
    echo "$meta.id" > ${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
