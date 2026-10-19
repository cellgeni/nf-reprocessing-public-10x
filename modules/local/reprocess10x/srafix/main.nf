process REPROCESS10X_SRAFIX {
    // Dumps an .sra that fastq-dump cannot, for runs REPROCESS10X_SRACHECK sent here.
    //
    // fastq-dump and fasterq-dump both need the virtual NAME column to write a record, and in
    // PRJNA524398's GenericFastq runs NAME is NULL for some spots ("row NULL while reading
    // column", error code 3, nothing written). The data columns are intact: vdb-validate
    // passes, and READ, QUALITY and READ_LEN read with vdb-dump. So the FASTQs are rebuilt from
    // those three: each spot's READ is cut by READ_LEN into one record per segment, written to
    // <run>_1/_2/_3.fastq.gz exactly as `fastq-dump --split-files` names them, with fastq-dump's
    // default defline (`@<run>.<spot> <spot> length=<len>`). Original names are lost; nothing
    // downstream reads them.
    //
    // vdb-dump is single-threaded, so the spot range is cut into one chunk per cpu, each
    // streamed through awk into its own gzip members; a file is the chunks' members in order,
    // which is a valid gzip stream. Nothing uncompressed reaches the disk.
    //
    // It refuses rather than emit something partial: every chunk must exit 0, every spot's
    // READ and QUALITY must be as long as its READ_LEN says, and every mate file must hold one
    // record per spot (a zero-length segment would drop records from one mate only).
    tag "${meta.id}"

    container "quay.io/cellgeni/reprocess_10x:latest"

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    spots=\$(vdb-dump --info ${sra} | awk '/^SEQ +:/{sub(/^SEQ +: */, ""); gsub(/,/, ""); print}')
    if ! [[ "\$spots" =~ ^[0-9]+\$ ]] || (( spots == 0 ))
    then
        echo "ERROR: could not read the spot count of ${prefix} (got '\$spots')." >&2
        exit 1
    fi
    chunks=${task.cpus}
    size=\$(( (spots + chunks - 1) / chunks ))
    echo "${prefix}: \$spots spots in \$chunks chunks of \$size"

    pids=()
    for (( c = 0; c < chunks; c++ ))
    do
        from=\$(( c * size + 1 ))
        to=\$(( (c + 1) * size ))
        if (( to > spots )); then to=\$spots; fi
        if (( from > to )); then break; fi
        (
            set -o pipefail
            vdb-dump -R "\$from-\$to" -f tab -C 'READ_LEN,READ,(INSDC:quality:text:phred_33)QUALITY' ${sra} \\
            | awk -F'\\t' -v acc="${prefix}" -v first="\$from" -v chunk="\$(printf '%04d' \$c)" '
                {
                    spot = first + NR - 1
                    n = split(\$1, len, /, */)
                    off = 1
                    for (i = 1; i <= n; i++) {
                        l = len[i] + 0
                        if (l > 0) {
                            out = "gzip -1 > part." chunk "_" i ".fastq.gz"
                            printf "@%s.%d %d length=%d\\n%s\\n+\\n%s\\n", acc, spot, spot, l, substr(\$2, off, l), substr(\$3, off, l) | out
                            count[i]++
                        }
                        off += l
                    }
                    if (length(\$2) != off - 1 || length(\$3) != off - 1) bad++
                    if (n > segs) segs = n
                }
                END {
                    for (i = 1; i <= segs; i++) {
                        close("gzip -1 > part." chunk "_" i ".fastq.gz")
                        print count[i] + 0 > ("part." chunk "_" i ".count")
                    }
                    print NR > ("part." chunk ".spots")
                    if (bad) {
                        printf "ERROR: %d spots in chunk %s whose READ or QUALITY length differs from READ_LEN\\n", bad, chunk > "/dev/stderr"
                        exit 1
                    }
                }'
        ) &
        pids+=( "\$!" )
    done
    status=0
    for pid in "\${pids[@]}"
    do
        wait "\$pid" || status=1
    done
    if (( status != 0 ))
    then
        echo "ERROR: a vdb-dump chunk of ${prefix} failed; refusing to emit a partial run." >&2
        exit 1
    fi

    dumped=\$(cat part.*.spots | awk '{s += \$1} END {print s + 0}')
    if (( dumped != spots ))
    then
        echo "ERROR: dumped \$dumped of \$spots spots of ${prefix}." >&2
        exit 1
    fi

    # one file per segment, chunks in order; every mate must hold one record per spot
    for seg in \$(ls part.*_*.count | sed -E 's/.*_([0-9]+)\\.count\$/\\1/' | sort -un)
    do
        records=\$(cat part.*_\${seg}.count | awk '{s += \$1} END {print s + 0}')
        echo "${prefix}_\${seg}.fastq.gz: \$records records"
        if (( records != spots ))
        then
            echo "ERROR: ${prefix}_\${seg} holds \$records records for \$spots spots: some spots have an empty segment." >&2
            echo "       Mates would disagree; refusing to emit a partial run." >&2
            exit 1
        fi
        cat \$(ls part.*_\${seg}.fastq.gz | sort) > "${prefix}_\${seg}.fastq.gz"
    done
    rm -f part.*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratoolkit: \$(vdb-dump --version | awk 'NF{print \$NF; exit}')
    END_VERSIONS
    """

    stub:
    """
    echo stub | gzip > ${meta.id}_1.fastq.gz
    echo stub | gzip > ${meta.id}_2.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratoolkit: stub
    END_VERSIONS
    """
}
