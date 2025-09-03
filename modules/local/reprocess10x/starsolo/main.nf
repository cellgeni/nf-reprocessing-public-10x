process REPROCESS10X_STARSOLO {
    tag "Renaming fastq files $prefix"

    container "docker://quay.io/cellgeni/reprocess_10x:latest"
    
    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    path(reference, stageAs: 'reference')
    path(whitelists, stageAs: 'whitelists')

    output:
    tuple val(meta), path("$prefix"), emit: mapping
    path "versions.yml", emit: versions

    script:
    bam_options = task.ext.bam_options ?: '--outSAMtype None --outReadsUnmapped Fastx'
    prefix = "${meta.id}"
    """
    # Move fastqs to sample directory
    mkdir -p "fastqs/${meta.id}"
    mv fastqs/*.gz "fastqs/${meta.id}/"

    # Run STARsolo
    source starsolo_10x_auto.sh
    starsolo_10x fastqs $prefix ${task.cpus} $reference $whitelists $bam_options

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version)
    END_VERSIONS
    """
}