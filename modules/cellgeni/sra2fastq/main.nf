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
    # fastq-dump drops a read entirely when the archive stores it with zero
    # length, writes only the surviving mate, and still exits 0. A single-mate
    # 10x run is useless, and left unchecked it travels two more stages before
    # surfacing as an unrelated-looking argument error in RENAME10XRUN. Keep the
    # dump log so the checks below can read what it actually did.
    parallel-fastq-dump -t $task.cpus -T . -s "${meta.id}" $args 2>&1 | tee dump.log

    if grep -qE 'Rejected [0-9]+ READS because READLEN < 1' dump.log
    then
        echo "ERROR: fastq-dump discarded a whole read of ${meta.id} for having zero length in the archive." >&2
        echo "       Only the surviving mate would have been emitted, which is not a usable 10x run." >&2
        grep -E 'Rejected [0-9]+ READS because READLEN < 1' dump.log | sort -u >&2
        exit 1
    fi

    shopt -s nullglob
    fastqs=( ${meta.id}*.fastq )
    if (( \${#fastqs[@]} < 2 ))
    then
        echo "ERROR: fastq-dump produced \${#fastqs[@]} FASTQ file(s) for ${meta.id}: \${fastqs[*]:-none}." >&2
        echo "       A 10x run needs at least a barcode read and a biological read." >&2
        exit 1
    fi

    # Every mate must hold the same number of records. Mismatched counts mean a
    # partial dump, and STARsolo would silently mispair reads.
    for i in "\${fastqs[@]}"
    do
        wc -l < "\$i" > "\$i.lines" &
    done
    wait

    expected=""
    for i in "\${fastqs[@]}"
    do
        lines=\$(cat "\$i.lines")
        # An empty count file means wc itself failed; arithmetic would read that
        # as zero and wave the file through.
        if [[ ! \$lines =~ ^[0-9]+\$ ]]
        then
            echo "ERROR: could not count records in \$i (wc reported '\$lines')." >&2
            exit 1
        fi
        if (( lines % 4 != 0 ))
        then
            echo "ERROR: \$i holds \$lines lines, which is not a whole number of FASTQ records." >&2
            exit 1
        fi
        records=\$(( lines / 4 ))
        echo "\$i: \$records records"
        if [[ -z \$expected ]]
        then
            expected=\$records
        elif [[ \$records != "\$expected" ]]
        then
            echo "ERROR: mates of ${meta.id} disagree on record count (\$i has \$records, expected \$expected)." >&2
            echo "       The dump is incomplete; refusing to emit a partial run." >&2
            exit 1
        fi
    done
    rm -f -- *.lines

    if (( expected == 0 ))
    then
        echo "ERROR: fastq-dump produced empty FASTQs for ${meta.id}." >&2
        exit 1
    fi

    # A bare 'wait' always returns 0, so a failed pigz used to pass unnoticed.
    # Wait on each PID individually and report the file that failed.
    pids=()
    for i in "\${fastqs[@]}"
    do
        pigz "\$i" &
        pids+=( "\$!" )
    done
    status=0
    for idx in "\${!pids[@]}"
    do
        if ! wait "\${pids[idx]}"
        then
            echo "ERROR: pigz failed to compress \${fastqs[idx]}" >&2
            status=1
        fi
    done
    if (( status != 0 )); then exit 1; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parallel-fastq-dump: \$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
        pigz: \$(pigz --version 2>&1 | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_1.fastq.gz
    touch ${meta.id}_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parallel-fastq-dump: \$(grep "fastq-dump" /versions.txt | cut -d ':' -f 2)
        pigz: \$(pigz --version 2>&1 | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
