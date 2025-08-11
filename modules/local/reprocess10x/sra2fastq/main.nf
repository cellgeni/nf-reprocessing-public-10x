process REPROCESS10X_SRA2FASTQ {
    tag "Converting sra file to .fastq for $meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x:latest':
        'quay.io/cellgeni/reprocess_10x:latest' }"

    input:
    tuple val(meta), path(sra)
    path whitelists

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"           , emit: versions

    script:
    """
    source sra_to_10x_fastq_gz.sh
    sra2fastq ${meta.id} ${whitelists} ${task.cpus}

    seqtk_version=\$(grep seqtk /versions.txt | cut -d ':' -f 2)
    fqdump_version=\$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk_version: \$seqtk_version
        parallel-fastq-dump: \$fqdump_version
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastq.gz

    seqtk_version=\$(grep seqtk /versions.txt | cut -d ':' -f 2)
    fqdump_version=\$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk_version: \$seqtk_version
        parallel-fastq-dump: \$fqdump_version
    END_VERSIONS
    """
}
