process REPROCESS10X_LOADDATA {
    tag "Loading $meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/reprocess_10x':
        'quay.io/cellgeni/reprocess_10x' }"

    input:
    tuple val(meta), val(links)

    output:
    tuple val(meta), path("${meta.id}.urls.list"), emit: urls
    tuple val(meta), path("SRR*"),                 emit: sra,   optional: true
    tuple val(meta), path(".f*q*"),                emit: fastq, optional: true
    tuple val(meta), path(".bam"),                 emit: bam,   optional: true

    script:
    def args = task.ext.args ?: ''
    def linkstring = links.join("\n")
    """
    # Create a file with links to download
    echo -e "$linkstring" > ${meta.id}.urls.list

    # Download data
    cat ${meta.id}.urls.list | xargs -I {} -n1 -P4 wget --retry-connrefused --read-timeout=20 --timeout=15 --tries=0 $args {}
    """

    stub:
    """
    touch ${meta.id}.urls.list
    touch stub.bam
    touch stub.fastq.gz
    touch SRRstub
    """
}
