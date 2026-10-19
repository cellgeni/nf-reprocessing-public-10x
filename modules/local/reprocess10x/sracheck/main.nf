process REPROCESS10X_SRACHECK {
    // Decides, per downloaded .sra, whether SRA2FASTQ can dump it or REPROCESS10X_SRAFIX has to.
    // Emits only a verdict, never the file: download10x joins the verdict back onto WGET10X's own
    // path, so SRA2FASTQ sees exactly the input it always did and its cached tasks stay cached.
    //
    // `fix` is returned for two signatures, both from PRJNA524398 (docs/archive-pathologies.md,
    // `null-name-column`):
    //   - fastq-dump stops with "row NULL while reading column": the virtual NAME column is NULL
    //     for some spots, and fastq-dump cannot write a record without it, whatever --defline.
    //     Probed on a few windows across the run rather than the whole file, which would cost as
    //     much as the dump itself. A NULL stretch the windows miss still fails in SRA2FASTQ, as
    //     every such run did before this step existed.
    //   - schema NCBI:SRA:GenericFastq: all 7 such runs seen so far either hit the NULL NAME or
    //     needed 65 GB to dump (SRR8816330). SRAFIX never reads NAME and streams in chunks.
    // Anything else, including a check that itself fails, is `ok`: the step may only ever add
    // a way through, never take one away.
    tag "${meta.id}"

    container "quay.io/cellgeni/reprocess_10x:latest"

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), env('VERDICT'), env('REASON'), emit: verdict
    path "versions.yml"                       , emit: versions

    script:
    def window = task.ext.window ?: 10000
    """
    VERDICT=ok
    REASON=""

    info=\$(vdb-dump --info ${sra} 2> info.err || true)
    # `SCHEMA : NCBI:SRA:GenericFastq:db#1` -- the value holds colons, so strip the label
    # rather than split on ':'
    spots=\$(echo "\$info" | awk '/^SEQ +:/{sub(/^SEQ +: */, ""); gsub(/,/, ""); print}')
    schema=\$(echo "\$info" | awk '/^SCHEMA +:/{sub(/^SCHEMA +: */, ""); print}')
    echo "spots=\${spots:-?} schema=\${schema:-?}"

    if [[ "\$schema" == NCBI:SRA:GenericFastq* ]]
    then
        VERDICT=fix
        REASON="schema \$schema"
    elif [[ "\$spots" =~ ^[0-9]+\$ ]] && (( spots > 0 ))
    then
        # windows at the start, the quartiles and the end
        for frac in 0 25 50 75 100
        do
            from=\$(( (spots - ${window}) * frac / 100 + 1 ))
            (( from < 1 )) && from=1
            to=\$(( from + ${window} - 1 ))
            (( to > spots )) && to=\$spots
            # the same options SRA2FASTQ dumps with, so a pass here is a pass there
            if ! fastq-dump -N "\$from" -X "\$to" -Z --split-files -F ${sra} > /dev/null 2> probe.err
            then
                if grep -q 'row NULL while reading column' probe.err
                then
                    VERDICT=fix
                    REASON="NULL NAME in spots \$from-\$to"
                    break
                fi
                echo "WARN: fastq-dump failed on spots \$from-\$to for another reason; leaving it to SRA2FASTQ:" >&2
                grep -v '^\$\\|=====' probe.err | head -3 >&2
            fi
        done
    else
        echo "WARN: could not read the spot count of ${meta.id}; leaving it to SRA2FASTQ" >&2
        cat info.err >&2 || true
    fi
    echo "verdict=\$VERDICT \$REASON"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratoolkit: \$(vdb-dump --version | awk 'NF{print \$NF; exit}')
    END_VERSIONS
    """

    stub:
    """
    VERDICT=ok
    REASON=stub
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratoolkit: stub
    END_VERSIONS
    """
}
