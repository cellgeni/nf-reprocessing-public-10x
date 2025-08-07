process REPROCESS10X_PARSEMETADATA {
    tag "Parsing metadata for $meta.id"
    container "${ workflow.containerEngine == 'singularity' ? 'docker://quay.io/cellgeni/reprocess_10x:dev': 'quay.io/cellgeni/reprocess_10x:dev' }"

    input:
    tuple val(meta), val(sample_ids)

    output:
    tuple val(meta), path("links.tsv"), emit: links
    tuple val(meta), path("*.list"), emit: list
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.txt"), emit: txt, optional: true
    tuple val(meta), path("*.soft"), emit: soft, optional: true
    path "versions.yml"           , emit: versions


    script:
    def datasetstring = sample_ids ? sample_ids.join("\n") : ''
    """
    # Create a sample list file
    echo -e "$datasetstring" > sample.list

    # Download and parse metadata
    collect_metadata.sh ${meta.id} ${sample_ids ? "sample.list" : ""}

    # Add sample IDs to metadata
    add_samples.awk ${meta.id}.sample_x_run.tsv ${meta.id}.parsed.tsv > links.tsv

    reprocess_version=\$(grep reprocess /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellgeni/reprocess_public_10x: \$reprocess_version
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    touch ${prefix}.bam

    reprocess_version=\$(grep reprocess /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellgeni/reprocess_public_10x: \$reprocess_version
    END_VERSIONS
    """
}
