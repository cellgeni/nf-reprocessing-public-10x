process REPROCESS10X_SRA2FASTQ {
    tag "Converting sra file to .fastq for $meta.id"
    cpus 16
    publishDir "results/$meta.dataset_id/raw/$meta.sample_id/", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x':
        'quay.io/cellgeni/reprocess_10x' }"

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

    reprocess_version=\$(grep reprocess /versions.txt | cut -d ':' -f 2)
    seqtk_version=\$(grep seqtk /versions.txt | cut -d ':' -f 2)
    fqdump_version=\$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellgeni/reprocess_public_10x: \$reprocess_version
        bamtofastq: \$bamtofastq_version
        parallel-fastq-dump: \$fqdump_version
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastq.gz

    reprocess_version=\$(grep reprocess /versions.txt | cut -d ':' -f 2)
    seqtk_version=\$(grep seqtk /versions.txt | cut -d ':' -f 2)
    fqdump_version=\$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellgeni/reprocess_public_10x: \$reprocess_version
        bamtofastq: \$bamtofastq_version
        parallel-fastq-dump: \$fqdump_version
    END_VERSIONS
    """
}
