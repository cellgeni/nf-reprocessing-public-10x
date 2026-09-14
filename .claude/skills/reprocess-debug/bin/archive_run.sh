#!/bin/bash
# Archive one reprocessing run off scratch, complete enough to read without
# Nextflow, a work dir, or access to Lustre.
#
#   archive_run.sh --batch 6
#   archive_run.sh --batch 5 --run tender_brattain
#   archive_run.sh --batch 6 --dry-run
#
# Writes into <dest>/batch<N>/<run_name>/:
#   failed<N>.log.gz       per-task stderr from the tools themselves — the file
#                          to read when you want the verbatim error
#   failures<N>.tsv        triage manifest (plain: small, and gets grepped)
#   failedjobs<N>.tsv      last failed attempt per task (plain)
#   triage<N>.json         triage.py --json aggregates
#   runlogs<N>.tsv.gz      full 36-field nextflow trace
#   batch<N>.csv           the input table, so the run is reproducible
#   execution_report_<ts>.html.gz, execution_timeline_<ts>.html.gz,
#   execution_trace_<ts>.txt.gz                     nextflow's own reports
#   lsf_{output,error}_<jobid>.log.gz               driver stdout/stderr
#   postmortem.html                                 the published write-up, if any;
#                          manifest.txt also carries its URL from docs/post-mortems.md
#   manifest.txt           what this run was, and what is and is not in here
#
# Nothing is skipped silently. A file that is missing or too large is recorded
# as such in manifest.txt, so a gap is visible from inside the archive.
set -euo pipefail

# /nfs/cellgeni/projects/ is owned by cellgeni-su at mode 755, so an ordinary
# group member cannot create under it. /nfs/cellgeni itself is group-writable.
DEST_ROOT="/nfs/cellgeni/reprocessing-runs"
BATCH=""; RUN=""; LATEST=0; DRY=0; FORCE=0; MAXSIZE=0   # MAXSIZE 0 = unlimited
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"

die()  { echo "error: $*" >&2; exit 1; }
warn() { echo "warn:  $*" >&2; }

usage() { sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    echo; echo "options: --batch N | --run NAME [--latest] [--dest DIR] [--dry-run]"
    echo "         [--max-size BYTES] [--force]"; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch)    BATCH="$2"; shift 2 ;;
        --run)      RUN="$2";   shift 2 ;;
        --dest)     DEST_ROOT="$2"; shift 2 ;;
        --max-size) MAXSIZE="$2"; shift 2 ;;
        --latest)   LATEST=1; shift ;;
        --dry-run)  DRY=1;    shift ;;
        --force)    FORCE=1;  shift ;;
        -h|--help)  usage ;;
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
STAMP="$(find_report_stamp "$(cut -f1 <<<"$RR_ROW")" 300 || true)"
JOBID="$(resolve_lsf_jobid "$RUN" || true)"
N="$BATCH"

echo "run   : $RUN (batch ${BATCH:-?})" >&2
awk -F'\t' '{printf "        started %s  status %s  ran %s  session %s\n", $1, $4, $2, $6}' <<<"$RR_ROW" >&2
echo "stamp : $STAMP    lsf job: ${JOBID:-not found}" >&2

DEST="$DEST_ROOT/batch${BATCH:-unknown}/$RUN"
echo "dest  : $DEST" >&2

# ------------------------------------------------------------------ what to copy
# The unsuffixed trio in data/ (failed.log, failedjobs.tsv, runlogs.tsv) is NOT
# batch 1's, despite living alongside it. failed.log and failedjobs.tsv are
# dated 2026-08-17 and the batch-1 run peaceful_feynman started 2026-08-19, so
# they describe an earlier run of the pre-`batches/` era and cannot be
# attributed from what is on disk. Archive them with archive_legacy.sh instead.
# Mapping --batch 1 onto them silently mislabels 1814 failures as this run's.
T_FAILED="data/tables/failed${N}.log"; T_JOBS="data/tables/failedjobs${N}.tsv"
T_LOGS="data/tables/runlogs${N}.tsv";  T_MANI="data/tables/failures${N}.tsv"

# src|dstname|gzip(0/1)
PLAN=()
add() { PLAN+=("$1|$2|$3"); }
add "$T_FAILED" "failed${N}.log"      1
add "$T_MANI"   "failures${N}.tsv"    0
add "$T_JOBS"   "failedjobs${N}.tsv"  0
add "$T_LOGS"   "runlogs${N}.tsv"     1
add "data/tables/batches/batch${N}.csv" "batch${N}.csv" 0
if [[ -n "$STAMP" ]]; then
    for kind in report timeline; do
        add "reports/execution_${kind}_${STAMP}.html" "execution_${kind}_${STAMP}.html" 1
    done
    add "reports/execution_trace_${STAMP}.txt" "execution_trace_${STAMP}.txt" 1
else
    NOTES_PRE="no reports/execution_report_* within 300s of the run's start time"
fi
if [[ -n "$JOBID" ]]; then
    add "logs/lsf/reprocessOutput${JOBID}.log" "lsf_output_${JOBID}.log" 1
    add "logs/lsf/reprocessError${JOBID}.log"  "lsf_error_${JOBID}.log"  1
fi
# the post-mortem source, whatever it was named this time
for pm in "data/failure-postmortem-batch${N}.html" "data/failure-postmortem-run${N}.html" \
          "data/run-${JOBID}-postmortem.html"; do
    [[ -f "$pm" ]] && { add "$pm" "postmortem.html" 0; break; }
done

# and the link to the published page, from the tracked index of record. Matched on
# the run name, which every batch row carries in its archive path. Absent is normal:
# the report is usually published after the run is archived.
PM_INDEX="$REPO/docs/post-mortems.md"
PM_URL=""
if [[ -f "$PM_INDEX" ]]; then
    PM_URL="$(grep -F "$RUN" "$PM_INDEX" 2>/dev/null \
        | grep -o 'https://claude\.ai/code/artifact/[0-9a-f-]\{36\}' | head -1 || true)"
fi

# triage json: generate next to the manifest if it is not there already
# --json takes the destination path; it does NOT write to stdout.
# Only batch 6 has a failures<N>.tsv — earlier runs were collected before the
# manifest format existed. triage.py reads those from the log pair instead, so
# an old run still gets its aggregates rather than being archived raw.
TRIAGE_SRC="data/tables/triage${N}.json"
if [[ ! -s "$TRIAGE_SRC" && $DRY -eq 0 ]]; then
    if [[ -s "$T_MANI" ]]; then
        echo "==> triage.py --manifest -> $TRIAGE_SRC" >&2
        "$BIN/triage.py" --manifest "$T_MANI" --json "$TRIAGE_SRC" >/dev/null \
            || warn "triage.py failed; archiving without triage${N}.json"
    elif [[ -s "$T_FAILED" && -s "$T_JOBS" ]]; then
        echo "==> triage.py --failed-log (no manifest for this run) -> $TRIAGE_SRC" >&2
        "$BIN/triage.py" --failed-log "$T_FAILED" --failedjobs "$T_JOBS" --json "$TRIAGE_SRC" >/dev/null \
            || warn "triage.py failed; archiving without triage${N}.json"
    fi
fi
# -s not -f: a failed run above leaves a zero-byte file, which is worse than none
[[ -s "$TRIAGE_SRC" ]] && add "$TRIAGE_SRC" "triage${N}.json" 0

# ------------------------------------------------------------------ report plan
hsize() { local b="${1:-0}"; awk -v b="$b" 'BEGIN{s="B KB MB GB TB";split(s,u," ");
    i=1; while(b>=1024 && i<5){b/=1024;i++} printf "%.1f%s", b, u[i]}'; }

# Destination name as the writer below will actually produce it. Keep this test
# identical to the one in the write loop: gzip only when the plan asks for it AND
# the file is over 1 MiB, so a small .log is archived plain, not as a .log.gz.
dst_name() { local dst="$1" gz="$2" sz="$3"
    if (( gz )) && (( sz > 1048576 )); then echo "$dst.gz"; else echo "$dst"; fi; }


# must be assigned, not just declared: `set -u` treats a declared-but-unset
# array as unbound when the plan turns out to have no gaps at all.
NOTES=()
[[ -n "${NOTES_PRE:-}" ]] && NOTES+=("$NOTES_PRE")
printf '\n%-52s %10s  %s\n' "SOURCE" "SIZE" "ACTION" >&2
for e in "${PLAN[@]}"; do
    IFS='|' read -r src dst gz <<<"$e"
    if [[ ! -f "$src" ]]; then
        printf '%-52s %10s  %s\n' "$src" "-" "MISSING — recorded in manifest" >&2
        NOTES+=("missing: $src (would have been $dst)"); continue
    fi
    sz=$(stat -c%s "$src")
    if (( MAXSIZE > 0 && sz > MAXSIZE )); then
        printf '%-52s %10s  %s\n' "$src" "$(hsize "$sz")" "SKIPPED — over --max-size" >&2
        NOTES+=("skipped: $src ($(hsize "$sz")) exceeded --max-size $(hsize "$MAXSIZE")")
        continue
    fi
    # Predict the destination name with the SAME test the writer uses below, or the
    # preview promises failed6.log.gz and the archive ends up holding failed6.log.
    printf '%-52s %10s  %s\n' "$src" "$(hsize "$sz")" "$(dst_name "$dst" "$gz" "$sz")" >&2
done

if (( DRY )); then
    echo >&2; echo "--dry-run: nothing written." >&2
    (( ${#NOTES[@]} )) && { echo "gaps:" >&2; printf '  %s\n' "${NOTES[@]}" >&2; }
    exit 0
fi

# ------------------------------------------------------------------ write
[[ -d "$DEST" && $FORCE -eq 0 ]] && die "$DEST exists — pass --force to overwrite"
mkdir -p "$DEST" || die "cannot create $DEST (is /nfs/cellgeni writable from here?)"

copied=0
for e in "${PLAN[@]}"; do
    IFS='|' read -r src dst gz <<<"$e"
    [[ -f "$src" ]] || continue
    sz=$(stat -c%s "$src")
    (( MAXSIZE > 0 && sz > MAXSIZE )) && continue
    if [[ "$(dst_name "$dst" "$gz" "$sz")" == *.gz ]]; then
        gzip -c "$src" > "$DEST/$dst.gz" && copied=$((copied+1))
    else
        cp "$src" "$DEST/$dst" && copied=$((copied+1))
    fi
done

# ------------------------------------------------------------------ manifest
{
    echo "run             : $RUN"
    echo "batch           : ${BATCH:-unknown}"
    awk -F'\t' '{printf "started         : %s\nduration        : %s\nnextflow status : %s\nsession         : %s\n", $1, $2, $4, $6}' <<<"$RR_ROW"
    echo "lsf driver job  : ${JOBID:-not found}"
    echo "report stamp    : $STAMP"
    echo "command         : $(cut -f7 <<<"$RR_ROW")"
    echo "nextflow        : $(nextflow -v 2>/dev/null || echo 'not on PATH at archive time')"
    echo "work dir root   : $REPO/nf-work"
    # the published page is the readable copy; postmortem.html here is the durable
    # one. Carry the link so the archive says where it is, not just that it exists.
    echo "post-mortem     : ${PM_URL:-not recorded in docs/post-mortems.md}"
    echo "archived        : $(date '+%Y-%m-%d %H:%M:%S') by ${USER:-unknown}"
    echo "archived from   : $REPO"
    echo
    echo "NOTE: nextflow status is not a verdict. errorStrategy 'ignore' keeps failed"
    echo "      tasks out of the exit status — a run recorded OK can have permanent"
    echo "      failures. Count them in failures<N>.tsv, and separate the ones that"
    echo "      self-healed on retry (column 'recovered')."
    echo
    echo "contents:"
    ( cd "$DEST" && ls -l --time-style=+ | awk 'NR>1{printf "  %-44s %s\n", $NF, $(NF-1)}' )
    if (( ${#NOTES[@]} )); then
        echo
        echo "not archived:"
        printf '  %s\n' "${NOTES[@]}"
    fi
} > "$DEST/manifest.txt"

echo >&2
echo "archived $copied files to $DEST ($(du -sh "$DEST" | cut -f1))" >&2
(( ${#NOTES[@]} )) && { echo "gaps recorded in manifest.txt:" >&2; printf '  %s\n' "${NOTES[@]}" >&2; }
echo "next: add the run to references/run-index.md" >&2
[[ -z "${PM_URL:-}" ]] && echo "      and, once a post-mortem is published, its URL to docs/post-mortems.md" >&2
