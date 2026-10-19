process REPROCESS10X_COLLATE {
    tag "run tables"

    // A biocontainer: Nextflow needs `ps` inside the image for task metrics, which python:*-slim lacks
    container "quay.io/biocontainers/python:3.12.12"

    input:
    path(datasets, stageAs: 'input/datasets.tsv')
    path(links, stageAs: 'links/links_*.tsv')
    path(sample_x_run, stageAs: 'sample_x_run/*')
    path(chemistry, stageAs: 'chemistry/*')
    path(soloqc, stageAs: 'soloqc/*')
    path(versions, stageAs: 'versions/versions_*.yml')
    val(samples)    // resolved sample metas, GroupKeys already unwrapped
    val(skipped)    // metas of samples dropped for want of a species

    output:
    path "tables/*",     emit: tables
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samples_json = groovy.json.JsonOutput.toJson(samples)
    def skipped_json = groovy.json.JsonOutput.toJson(skipped)
    """
    cat <<'END_JSON' > samples.json
    ${samples_json}
    END_JSON
    cat <<'END_JSON' > skipped.json
    ${skipped_json}
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d ' ' -f 2)
    END_VERSIONS

    collate_run_tables.py \\
        --datasets input/datasets.tsv \\
        --links links \\
        --sample-x-run sample_x_run \\
        --chemistry chemistry \\
        --soloqc soloqc \\
        --versions versions \\
        --extra-versions versions.yml \\
        --samples samples.json \\
        --skipped skipped.json \\
        --out tables
    """

    stub:
    """
    mkdir -p tables
    touch tables/versions.yml tables/datasets.tsv tables/samples.tsv tables/runs.tsv tables/run_chemistry.tsv tables/mapping_qc_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
