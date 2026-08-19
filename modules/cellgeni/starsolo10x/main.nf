/*
 * Module: cellgeni/starsolo
 */

process STARSOLO10X {
    tag "${meta.id}"
    container "docker://quay.io/cellgeni/starsolo:v4.2"
    
    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    tuple val(ref_meta), path(reference, stageAs: "reference")

    output:
    tuple val(meta), path("${meta.id}"), emit: mapping
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def args   = task.ext.args ?: ""
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
    starsolo 10x fastqs "$prefix" --ref "\$PWD/${ref_meta.id}/reference" --cpus ${task.cpus} $args

    cat <<-END_VERSIONS > "\$workdir/versions.yml"
    "${task.process}":
        STAR: \$(STAR --version)
        cellgeni/STARsolo: \$(starsolo --version)
    END_VERSIONS
    """
}
