def getFileName(filepath, sample_id) {
    if ( filepath ==~ /.*fastq\.gz/ ) {
        return "fastq/${sample_id}/$filepath"
    } else if ( filepath ==~ /.*\.bam/ ) {
        return "bam/${sample_id}/$filepath"
    } else if ( filepath != "versions.yml" ) {
        return "sra/${sample_id}/$filepath"
    } else {
        return null
    }
}

process REPROCESS10X_LOADDATA {
    tag "Loading ${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x:latest':
        'quay.io/cellgeni/reprocess_10x:latest' }"

    input:
    tuple val(meta), val(link)

    output:
    tuple val(meta), path("${meta.id}"),             emit: sra,   optional: true
    tuple val(meta), path("*.f*q*"),               emit: fastq, optional: true
    tuple val(meta), path("*.bam"),                emit: bam,   optional: true
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"
    """
    # Download file
    wget $args $link

    # Rename downloaded files
    shopt -s extglob
    if [[ "${meta.type}" == "BAM" && ! -f "${prefix}.bam" ]]; then
        mv -T *.bam* "${prefix}.bam"
    elif [[ "${meta.type}" == "SRA" && ! -f "${prefix}" ]]; then
        mv -T SRR!(*.urls.list) "${prefix}"
    fi

    # save versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | awk '{ print \$3 }')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.urls.list
    touch ${prefix}.bam
    touch ${prefix}_1.fastq.gz
    touch ${prefix}_2.fastq.gz
    touch ${prefix}

    # save versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | awk '{ print \$3 }')
    END_VERSIONS
    """
}
