process WGET10X {
    // One task fetches one URL, and a run whose sample sits in several datasets is
    // fetched once per dataset, so "${meta.id}" alone is not unique: the Aug 2026
    // batch had four concurrent tasks sharing the tag "Loading SRR20167168". Tags
    // are the key every triage script joins on, so this one carries the file name
    // and avoids spaces — a space is what made the block-splitter regex in
    // docs/agent_debug.md drop every WGET10X failure it was asked to read.
    tag "${meta.id}:${(link?.toString()?.trim()?.tokenize('/') ?: ['no-url'])[-1]}"

    container "quay.io/cellgeni/reprocess_10x:latest"

    input:
    tuple val(meta), val(link)

    output:
    tuple val(meta), path("${meta.id}"),             emit: sra,   optional: true
    tuple val(meta), path("*.f*q*"),               emit: fastq, optional: true
    tuple val(meta), path("*.bam"),                emit: bam,   optional: true
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def verify = task.ext.verify_downloads == null ? true : task.ext.verify_downloads
    def tries  = task.ext.dl_tries ?: 6
    """
    # An empty URL reaches wget as no argument at all, which makes it print its
    # usage and exit 1 — an opaque failure for what is really a metadata gap.
    if [[ -z "\$(echo '$link' | tr -d '[:space:]')" ]]
    then
        echo "ERROR: no download URL for ${meta.id} (type ${meta.type})." >&2
        echo "       The metadata row for this run has an empty location field." >&2
        exit 1
    fi

    # Stagger Nextflow-level retries so they do not land in the same throttling
    # window the previous attempt was rejected in. wget handles backoff within an
    # attempt; this spaces the attempts themselves apart.
    if (( ${task.attempt} > 1 ))
    then
        backoff=\$(( 30 * (${task.attempt} - 1) * (${task.attempt} - 1) ))
        echo "Attempt ${task.attempt}: waiting \${backoff}s before retrying $link" >&2
        sleep "\$backoff"
    fi

    # Download, resuming this task's own partial between in-job attempts.
    #
    # ENA drops long transfers mid-stream: in the Aug 2026 batch every download
    # failure was a "Read error … (Connection timed out)", and the two runs that
    # were lost outright were 33 GB files that never finished. A Nextflow-level
    # retry cannot rescue those, because each attempt starts in a clean work dir
    # and restarts from byte zero — five attempts at a file that dies after a few
    # gigabytes make no more progress than one. Resuming in-job is what turns an
    # unreliable link into a completed file.
    #
    # --continue therefore resumes only the partial this task wrote itself, never
    # one inherited from elsewhere, and the integrity check below is what makes
    # that safe: a resumed file that spliced badly is caught and refetched clean
    # rather than handed downstream.
    rc=0
    for attempt in \$(seq 1 ${tries})
    do
        rc=0
        wget $args --continue $link || rc=\$?
        if (( rc == 0 ))
        then
            break
        fi
        if (( attempt < ${tries} ))
        then
            pause=\$(( 20 * attempt ))
            echo "NOTE: wget exited \$rc on in-job attempt \$attempt/${tries} for $link; resuming in \${pause}s" >&2
            sleep "\$pause"
        fi
    done
    if (( rc != 0 ))
    then
        echo "ERROR: download of $link did not complete after ${tries} in-job attempts (last wget exit \$rc)." >&2
        echo "       Partial data is left in the work dir; the run was not passed downstream." >&2
        exit 4
    fi

    # One task downloads one URL, so exactly one file must have appeared.
    shopt -s nullglob
    mapfile -t downloaded < <(find . -maxdepth 1 -type f ! -name '.*' ! -name 'versions.yml' -printf '%f\\n')
    if (( \${#downloaded[@]} != 1 ))
    then
        echo "ERROR: expected exactly 1 file from $link, found \${#downloaded[@]}: \${downloaded[*]:-none}" >&2
        exit 1
    fi
    got=\${downloaded[0]}
    if [[ ! -s "\$got" ]]
    then
        echo "ERROR: \$got downloaded from $link is empty" >&2
        exit 1
    fi

    ${verify ? """
    # Verify the archive decompresses cleanly. wget accepts a truncated body when
    # the server sends no Content-Length, and it also reports success for a file
    # that reached its advertised byte count but spliced wrongly across a resume —
    # 11 downloads in the Aug 2026 batch arrived at exactly the expected size and
    # still failed CRC, so a length comparison would have passed every one. Left
    # undetected these surface stages later as gzip CRC or end-of-stream errors.
    archive_ok() {
        local f="\$1"
        case "\$f" in
            *.gz|*.bgz)
                pigz -t -p ${task.cpus} "\$f" 2>gzip_test.err || gzip -t "\$f" 2>>gzip_test.err
                ;;
            *.bz2)
                bzip2 -t "\$f" 2>gzip_test.err
                ;;
            *)
                return 0
                ;;
        esac
    }

    if ! archive_ok "\$got"
    then
        # A bad archive after a resumed transfer is most likely a bad splice, so
        # discard it and fetch the whole file once, cleanly. Only a second failure
        # means the source itself is serving corrupt data.
        echo "NOTE: \$got failed its integrity check; discarding it and refetching $link from scratch." >&2
        cat gzip_test.err >&2
        rm -f "\$got" gzip_test.err
        wget $args $link
        if ! archive_ok "\$got"
        then
            echo "ERROR: \$got is not a valid archive; the download from $link is corrupt or truncated." >&2
            echo "       This failed on a clean full download as well as a resumed one." >&2
            cat gzip_test.err >&2
            exit 1
        fi
    fi
    rm -f gzip_test.err
    """ : ''}

    # Rename downloaded files
    shopt -s extglob
    if [[ "${meta.type}" == "BAM" && ! -f "${prefix}.bam" ]]; then
        mv -T *.bam* "${prefix}.bam"
    elif [[ "${meta.type}" == "SRA" && ! -f "${prefix}" ]]; then
        mv -T [ESD]RR!(*.urls.list) "${prefix}"
    fi

    # save versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | awk '{ print \$3 }')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    # One task downloads exactly one URL. Creating a file per output type instead
    # made every stub task emit on the sra, fastq and bam channels at once, which
    # does not happen in a real run and left the downstream channel shapes
    # untestable. Derive the single file from the link, as the real script does.
    target=\$(basename "$link")
    touch "\$target"

    shopt -s extglob
    if [[ "${meta.type}" == "BAM" && ! -f "${prefix}.bam" ]]; then
        mv -T *.bam* "${prefix}.bam"
    elif [[ "${meta.type}" == "SRA" && ! -f "${prefix}" ]]; then
        mv -T [ESD]RR!(*.urls.list) "${prefix}"
    fi

    # save versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | awk '{ print \$3 }')
    END_VERSIONS
    """
}
