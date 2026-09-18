#!/bin/bash
# Copy a run's LSF driver stdout/stderr onto NFS right away, without waiting
# for the full archive_run.sh pass.
#
#   save_lsf_log.sh --batch 7
#   save_lsf_log.sh --run elegant_lamarr
#
# archive_run.sh already copies these two files (§8 in SKILL.md), but only as
# part of archiving the whole run — which happens "when you are done" with a
# batch, sometimes much later, sometimes never. logs/lsf/ lives on scratch,
# which is not backed up, so the driver log is one `lfs migrate`/cleanup away
# from being unrecoverable in the meantime. This script does just the log copy,
# safe to run immediately after a debug session and safe to re-run — it never
# refuses because the destination exists, unlike archive_run.sh's full pass,
# and archiving the run later just overwrites these two files with the same
# bytes.
#
# Writes into <dest>/batch<N>/<run_name>/ (same tree archive_run.sh uses):
#   lsf_output_<jobid>.log[.gz]   lsf_error_<jobid>.log[.gz]
set -euo pipefail

DEST_ROOT="/nfs/cellgeni/reprocessing-runs"
BATCH=""; RUN=""; LATEST=0
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"

die() { echo "error: $*" >&2; exit 1; }

usage() { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    echo; echo "options: --batch N | --run NAME [--latest] [--dest DIR]"; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch)  BATCH="$2"; shift 2 ;;
        --run)    RUN="$2";   shift 2 ;;
        --dest)   DEST_ROOT="$2"; shift 2 ;;
        --latest) LATEST=1; shift ;;
        -h|--help) usage ;;
        *) die "unknown argument: $1 (try --help)" ;;
    esac
done
[[ -n "$BATCH" || -n "$RUN" ]] || usage

cd "$REPO"
[[ -f main.nf ]] || die "no main.nf in $REPO — wrong repo root?"
BIN="$(dirname "${BASH_SOURCE[0]}")"
# shellcheck source=resolve_run.sh
source "$BIN/resolve_run.sh"

resolve_run "$BATCH" "$RUN" "$LATEST" || exit 1
RUN="$RR_RUN"; BATCH="$RR_BATCH"
JOBID="$(resolve_lsf_jobid "$RUN" || true)"
[[ -n "$JOBID" ]] || die "no logs/lsf/reprocessOutput*.log mentions [$RUN] — nothing to save"

DEST="$DEST_ROOT/batch${BATCH:-unknown}/$RUN"
mkdir -p "$DEST" || die "cannot create $DEST (is /nfs/cellgeni writable from here?)"

copied=0
for pair in "reprocessOutput${JOBID}.log:lsf_output_${JOBID}.log" \
            "reprocessError${JOBID}.log:lsf_error_${JOBID}.log"; do
    src="logs/lsf/${pair%%:*}"; dst="${pair##*:}"
    if [[ ! -f "$src" ]]; then
        echo "warn:  $src not found, skipping" >&2
        continue
    fi
    sz=$(stat -c%s "$src")
    if (( sz > 1048576 )); then
        gzip -c "$src" > "$DEST/$dst.gz" && { copied=$((copied+1)); echo "saved  $DEST/$dst.gz" >&2; }
    else
        cp "$src" "$DEST/$dst" && { copied=$((copied+1)); echo "saved  $DEST/$dst" >&2; }
    fi
done

(( copied > 0 )) || die "found jobid $JOBID but neither log file exists in logs/lsf/"
echo "saved $copied LSF log file(s) for $RUN to $DEST" >&2
