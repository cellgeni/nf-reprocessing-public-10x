process REPROCESS10X_LOADDATA {
    tag "Loading $meta.id"
    publishDir "results/$meta.dataset_id/raw/$meta.sample_id/", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x':
        'quay.io/cellgeni/reprocess_10x' }"

    input:
    tuple val(meta), val(links)

    output:
    tuple val(meta), path("${meta.id}.urls.list"), emit: urls
    tuple val(meta), path("$meta.id"),             emit: sra,   optional: true
    tuple val(meta), path("*.f*q*"),               emit: fastq, optional: true
    tuple val(meta), path("*.bam"),                emit: bam,   optional: true

    when:

    script:
    def args = task.ext.args ?: ''
    def linkstring = links.join("\n")
    def prefix = "${meta.id}"
    """
    # Create a file with links to download
    echo -e "$linkstring" > ${meta.id}.urls.list

    # Download data
    cat ${meta.id}.urls.list | xargs -I {} -n1 -P4 wget --retry-connrefused --read-timeout=20 --timeout=15 --tries=0 $args {}

    # Rename downloaded files
    if "${meta.type}" == "BAM"; then
        mv -T *.bam* "${prefix}.bam"
    elif "${meta.type}" == "SRA"; then
        mv -T SRR* "${prefix}"
    fi
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.urls.list
    touch ${prefix}.bam
    touch ${prefix}_1.fastq.gz
    touch ${prefix}_2.fastq.gz
    touch ${prefix}
    """
}
