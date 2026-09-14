/*
 * Module: cellgeni/rename10xrun
 */

process RENAME10XRUN {
    tag "${meta.id}"
    container "quay.io/cellgeni/metacells-python:latest"

    input:
    tuple val(meta), path(fastqs, stageAs: "fastqs/*")
    path  whitelist_dir

    output:
    // Index reads travel with the biological ones. `reads` is what the run is
    // carried by downstream — it is joined with the chemistry report and handed to
    // RENAME10XSAMPLE — so an I1/I2 file that only ever appeared on the optional
    // `index` channel was dropped outright: nothing consumed that channel, and the
    // run's index reads were missing from both the renamed sample set and the
    // published per-run output, which the README documents as holding I1.
    // RENAME10XSAMPLE is already built for them (its ROLE_RE matches I1/I2, it
    // handles multi-part _I1_002 files, and it treats index length mismatches as
    // non-fatal); they simply never arrived. `index` is kept as a narrower view for
    // any consumer that wants only the index reads.
    // Note the sample-level split still holds: RENAME10XSAMPLE emits R files on its
    // own `reads` and index files on its own `index`, so the aligners are unaffected.
    tuple val(meta), path("*_[RI]*_001.fastq.gz"), emit: reads
    tuple val(meta), path("*_I*_001.fastq.gz"), optional: true, emit: index
    tuple val(meta), path("${meta.id}.chemistry.json"), emit: chemistry
    path  "versions.yml",                               emit: versions

    script:
    def args = task.ext.args ?: ""
    """
    infer_10x_run_recommended.py \\
        --fastqs ${fastqs} \\
        --run-id ${meta.id} \\
        --whitelist-dir ${whitelist_dir} \\
        ${args} \\
        --outdir . \\
        --json ${meta.id}.chemistry.json \\
        --tsv ${meta.id}.chemistry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def lane = meta.lane ?: 1
    """
    echo "stub" | gzip > ${meta.id}_S1_L1_R1_001.fastq.gz
    echo "stub" | gzip > ${meta.id}_S1_L1_R2_001.fastq.gz
    echo '{"sample":"${meta.id}","chemistry":"stub"}' > ${meta.id}.chemistry.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
