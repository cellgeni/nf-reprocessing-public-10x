process SRA2FASTQ {
    tag "${meta.id}"

    container "quay.io/cellgeni/reprocess_10x:latest"

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: '--split-files -F'
    """
    parallel-fastq-dump -t $task.cpus -T . -s "${meta.id}" $args

    for i in ${meta.id}*fastq
    do
        [[ -e "\$i" ]] || continue
        pigz "\$i" &
    done
    wait

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parallel-fastq-dump: \$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
        pigz: \$(pigz --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_1.fastq.gz
    touch ${meta.id}_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parallel-fastq-dump: \$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
        pigz: \$(pigz --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
