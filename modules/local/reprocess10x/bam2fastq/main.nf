process REPROCESS10X_BAM2FASTQ {
    tag "Converting .bam file to .fastq for $meta.id"
    cpus 16
    publishDir "results/$meta.dataset_id/raw/$meta.sample_id/", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x:latest':
        'quay.io/cellgeni/reprocess_10x:latest' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"           , emit: versions

    script:
    """
    source bam_to_10x_fastq_gz.sh
    bam2fastq ${meta.id} ${task.cpus}

    reprocess_version=\$(grep reprocess /versions.txt | cut -d ':' -f 2)
    bamtofastq_version=\$(grep bamtofastq /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellgeni/reprocess_public_10x: \$reprocess_version
        bamtofastq: \$bamtofastq_version
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastq.gz

    reprocess_version=\$(grep reprocess /versions.txt | cut -d ':' -f 2)
    bamtofastq_version=\$(grep bamtofastq /versions.txt | cut -d ':' -f 2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellgeni/reprocess_public_10x: \$reprocess_version
        bamtofastq: \$bamtofastq_version
    END_VERSIONS
    """
}
