process REPROCESS10X_RUNREPORT {
    tag "run report"
    label 'process_single'

    // MultiQC ships its own python3, which the stdlib-only run_report.py runs on
    container "quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2"

    input:
    path(tables, stageAs: 'tables/*')       // REPROCESS10X_COLLATE's run-level tables
    path(starsolo, stageAs: 'starsolo/*')   // STARsolo sample dirs, named by sample id
    path(cellranger, stageAs: 'cellranger/*') // Cell Ranger outputs: not parsed yet, they make this wait for them
    path(trace)                             // this run's own trace, complete but for this task
    val(run_info)                           // workflow runName, sessionId, commitId, commandLine, start
    val(predicted)                          // [sample, reads, mem_gb, hours] from meta.resources

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def run_info_json  = groovy.json.JsonOutput.toJson(run_info)
    def predicted_json = groovy.json.JsonOutput.toJson(predicted)
    """
    cat <<'END_JSON' > run_info.json
    ${run_info_json}
    END_JSON
    cat <<'END_JSON' > predicted.json
    ${predicted_json}
    END_JSON

    run_report.py \\
        --trace ${trace} \\
        --tables tables \\
        --starsolo starsolo \\
        --run-info run_info.json \\
        --predicted predicted.json \\
        --out mqc

    # MultiQC builds multiqc_data/ under TMPDIR and copies it across with its modes, which
    # drops the work dir's setgid bit: files then land in the user's primary group, which
    # has no quota on /lustre ("Disk quota exceeded"). A TMPDIR inside the work dir keeps
    # everything in the inherited group.
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"
    multiqc --force --config mqc/multiqc_config.yaml mqc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(python3 -c "import importlib.metadata as m; print(m.version('multiqc'))")
        python: \$(python3 --version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(python3 -c "import importlib.metadata as m; print(m.version('multiqc'))")
        python: \$(python3 --version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
