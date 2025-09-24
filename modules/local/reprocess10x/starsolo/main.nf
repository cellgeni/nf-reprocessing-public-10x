process REPROCESS10X_STARSOLO {
    tag "Renaming fastq files $prefix"

    container "docker://quay.io/cellgeni/reprocess_10x:latest"
    
    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    tuple val(ref_meta), path(reference, stageAs: "reference")
    path(whitelists, stageAs: 'whitelists')

    output:
    tuple val(meta), path("$prefix"), emit: mapping
    path "versions.yml", emit: versions

    script:
    bam_options = task.ext.bam_options ?: '--outSAMtype None --outReadsUnmapped Fastx'
    prefix = "${meta.id}"
    """
    # Write the path to workdir
    workdir=\$PWD

    # Move fastqs to sample directory
    mkdir -p "fastqs/${meta.id}"
    mv fastqs/*.gz "fastqs/${meta.id}/"

    # Rename reference directory so that it contains specie name
    mkdir -p "${ref_meta.id}"
    mv "$reference" "${ref_meta.id}/reference"    

    # Run STARsolo
    source starsolo_10x_auto.sh
    starsolo_10x fastqs $prefix ${task.cpus} "\$PWD/${ref_meta.id}/reference" "\$PWD/$whitelists" "$bam_options" "${ref_meta.id}"

    cat <<-END_VERSIONS > "\$workdir/versions.yml"
    "${task.process}":
        STAR: \$(STAR --version)
    END_VERSIONS
    """
}