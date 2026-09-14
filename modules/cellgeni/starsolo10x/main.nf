/*
 * Module: cellgeni/starsolo
 */

process STARSOLO10X {
    tag "${meta.id}"
    // v4.3 is the first wrapper release with `--wl` and `--skip-length-checks`.
    container "docker://quay.io/cellgeni/starsolo:v4.3"
    
    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    tuple val(ref_meta), path(reference, stageAs: "reference")

    output:
    tuple val(meta), path("${meta.id}"), emit: mapping
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def args   = task.ext.args ?: ""
    // meta.wl holds a chemistry id from upstream run-level inference, in the very
    // vocabulary `starsolo 10x --wl` speaks. Given it, the wrapper takes the whitelist
    // and CB/UMI geometry from the id instead of matching barcodes against every
    // whitelist in a subsample — which is where runs of perfectly good public data were
    // being lost, since barcode matching reads the head of the FASTQ.
    //
    // --skip-length-checks travels with it, and only with it: the read lengths behind
    // that id were already validated per run (variable-length barcode reads, minimum
    // biological read length, CB/UMI geometry) and again per sample when the runs were
    // merged, so re-deciding them here from a fresh head sample only adds a way to
    // disagree. It downgrades those aborts to warnings; the wrapper still trims the UMI
    // to fit R1, and still fails outright if R1 cannot hold the barcode. A sample with
    // no inference to stand on (BAM-derived runs carry no report) keeps every check.
    //
    // Both flags precede ext.args so that a `--` passthrough there still terminates the
    // wrapper's own options, and an explicit --wl in a config still wins (last wins).
    def chemistry = meta.wl ? "--wl ${meta.wl} --skip-length-checks" : ""
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
    starsolo 10x fastqs "$prefix" --ref "\$PWD/${ref_meta.id}/reference" --cpus ${task.cpus} $chemistry $args

    cat <<-END_VERSIONS > "\$workdir/versions.yml"
    "${task.process}":
        STAR: \$(STAR --version)
        cellgeni/STARsolo: \$(starsolo --version)
    END_VERSIONS
    """
}
