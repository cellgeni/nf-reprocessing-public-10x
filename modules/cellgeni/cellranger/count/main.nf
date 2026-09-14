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
    // meta.chemistry is set upstream only when run-level inference identified a chemistry
    // that Cell Ranger's own --chemistry=auto will not accept on its own (Multiome
    // ARC-v1). It replaces rather than appends to any --chemistry already in ext.args,
    // since Cell Ranger takes the flag once; when it is unset, args is passed through
    // untouched, so an explicit --chemistry in a config still wins for every other sample.
    if (meta.chemistry) {
        args = (args.replaceAll(/--chemistry[=\s]+\S+/, '').trim() + " --chemistry=${meta.chemistry}").trim()
    }
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
