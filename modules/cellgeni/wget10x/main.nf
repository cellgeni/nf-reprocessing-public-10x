process WGET10X {
    tag "Loading ${meta.id}"

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

    # Download file
    wget $args $link

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
    # the server sends no Content-Length, and the resulting partial gzip is only
    # noticed further down the pipeline as a CRC or end-of-stream error.
    case "\$got" in
        *.gz|*.bgz)
            if ! pigz -t -p ${task.cpus} "\$got" 2>gzip_test.err && ! gzip -t "\$got" 2>>gzip_test.err
            then
                echo "ERROR: \$got is not a valid gzip stream; the download from $link is corrupt or truncated." >&2
                cat gzip_test.err >&2
                exit 1
            fi
            rm -f gzip_test.err
            ;;
        *.bz2)
            if ! bzip2 -t "\$got"
            then
                echo "ERROR: \$got is not a valid bzip2 stream; the download from $link is corrupt or truncated." >&2
                exit 1
            fi
            ;;
    esac
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
