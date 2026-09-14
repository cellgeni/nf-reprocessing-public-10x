#!/bin/bash
# Collect the failure artifacts for one reprocessing run.
#
#   collect_run_logs.sh --batch 6
#   collect_run_logs.sh --run spontaneous_ampere --suffix 6
#
# Writes, into --outdir (default data/tables):
#   runlogs<S>.tsv     full 36-field nextflow trace export, header included
#   failedjobs<S>.tsv  last FAILED attempt per task: attempt exit name tag hash workdir
#   failed<S>.log      per-task stderr in "### <tag> (<workdir>) ###" blocks,
#                      LSF report stripped, TERM_* kill reason kept as one line
#   failures<S>.tsv    manifest: process tag exit attempts_failed recovered workdir first_error
#
# failures<S>.tsv is what bin/triage.py wants. The other three exist because
# earlier runs were collected by hand and the formats are load-bearing.
set -euo pipefail

FIELDS="accelerator,accelerator_type,attempt,complete,container,cpu_model,cpus,disk,duration,exit,hash,hostname,memory,module,name,native_id,pcpu,peak_rss,peak_vmem,pmem,process,queue,read_bytes,realtime,rss,scratch,start,status,submit,tag,task_id,time,vmem,wchar,workdir,write_bytes"

BATCH=""; RUN=""; SUFFIX=""; OUTDIR="data/tables"; FORCE=0; LATEST=0
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"

die()  { echo "error: $*" >&2; exit 1; }
warn() { echo "warn:  $*" >&2; }

usage() {
    sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    echo
    echo "options: --batch N | --run NAME [--suffix S] [--outdir DIR] [--latest] [--force]"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch)  BATCH="$2"; shift 2 ;;
        --run)    RUN="$2";   shift 2 ;;
        --suffix) SUFFIX="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --latest) LATEST=1; shift ;;
        --force)  FORCE=1;  shift ;;
        -h|--help) usage ;;
        *) die "unknown argument: $1 (try --help)" ;;
    esac
done
[[ -n "$BATCH" || -n "$RUN" ]] || usage

cd "$REPO"
[[ -f main.nf && -f .nextflow/history ]] || die "no .nextflow/history in $REPO — wrong repo root?"

# ---------------------------------------------------------------- tooling
# nextflow log reads .nextflow/cache, written by the version that ran the
# pipeline. The PATH default on the farm is 25.04.4 and cannot read a 26.x
# cache; the repo manifest requires >=26.04.1 anyway.
load_module() {
    [[ -z "${MODULESHOME:-}" && -f /etc/profile.d/modules.sh ]] && . /etc/profile.d/modules.sh
    if command -v module >/dev/null 2>&1 || [[ "$(type -t module || true)" == function ]]; then
        module load "$1" >/dev/null 2>&1 || true
    fi
}
ensure_tool() {
    local bin="$1" mod="$2"
    command -v "$bin" >/dev/null 2>&1 && return 0
    load_module "$mod"
    command -v "$bin" >/dev/null 2>&1 || die "$bin not on PATH. Run: module load $mod"
}
ensure_tool duckdb   cellgen/duckdb
ensure_tool nextflow cellgen/nextflow/26.04.6

# The farm PATH default is 25.04.4, which reads a 26.x .nextflow/cache as an
# empty run — no error, no rows. Having *a* nextflow is not enough.
nf_ver() { nextflow -v 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true; }
if [[ "$(nf_ver)" != 26.* ]]; then
    warn "nextflow $(nf_ver) cannot read a 26.x cache — loading cellgen/nextflow/26.04.6"
    load_module cellgen/nextflow/26.04.6
    hash -r 2>/dev/null || true
fi
[[ "$(nf_ver)" == 26.* ]] || die "need nextflow 26.x, have $(nf_ver). Run: module load cellgen/nextflow/26.04.6"

# ---------------------------------------------------------------- resolve run
# .nextflow/history: timestamp \t duration \t run_name \t status \t script_id \t session \t command
if [[ -z "$RUN" ]]; then
    mapfile -t cands < <(awk -F'\t' -v b="batch${BATCH}.csv" 'index($7, b) {print}' .nextflow/history)
    (( ${#cands[@]} )) || die "no run in .nextflow/history mentions batch${BATCH}.csv"
    if (( ${#cands[@]} > 1 && LATEST == 0 )); then
        echo "batch $BATCH has ${#cands[@]} runs (a -resume keeps the session and takes a new name):" >&2
        printf '%s\n' "${cands[@]}" | awk -F'\t' '{printf "  %s  %-22s %-4s %s\n", $1, $3, $4, $2}' >&2
        die "pick one with --run <name>, or take the most recent with --latest"
    fi
    row="${cands[-1]}"
    RUN="$(cut -f3 <<<"$row")"
    echo "resolved batch $BATCH -> $RUN" >&2
    awk -F'\t' '{printf "  started %s  status %s  ran %s  session %s\n", $1, $4, $2, $6}' <<<"$row" >&2
else
    row="$(awk -F'\t' -v r="$RUN" '$3 == r {print}' .nextflow/history | tail -1)"
    [[ -n "$row" ]] || die "run '$RUN' not in .nextflow/history"
    [[ -z "$BATCH" ]] && BATCH="$(grep -oE 'batch[0-9]+\.csv' <<<"$row" | head -1 | grep -oE '[0-9]+' || true)"
fi
[[ -n "$SUFFIX" ]] || SUFFIX="${BATCH:-$RUN}"

mkdir -p "$OUTDIR"
TRACE="$OUTDIR/runlogs${SUFFIX}.tsv"
JOBS="$OUTDIR/failedjobs${SUFFIX}.tsv"
BLOB="$OUTDIR/failed${SUFFIX}.log"
MANI="$OUTDIR/failures${SUFFIX}.tsv"
for f in "$TRACE" "$JOBS" "$BLOB" "$MANI"; do
    [[ -e "$f" && $FORCE -eq 0 ]] && die "$f exists — pass --force to overwrite"
done

TMPD="$(mktemp -d)"; trap 'rm -rf "$TMPD"' EXIT

# ---------------------------------------------------------------- 1. trace
echo "==> nextflow log $RUN -> $TRACE" >&2
tr ',' '\t' <<<"$FIELDS" > "$TRACE"
nextflow log "$RUN" -f "$FIELDS" >> "$TRACE"
echo "    $(( $(wc -l < "$TRACE") - 1 )) task rows" >&2

# ---------------------------------------------------------------- 2. failed jobs
# One row per task name: its last FAILED attempt, plus how many attempts failed
# and whether any attempt of that task later succeeded. Do NOT read `attempt`
# as diagnostic — nextflow.config's errorStrategy retries everything once.
cat > "$TMPD/q.sql" <<SQL
CREATE TEMP TABLE t AS
  SELECT * FROM read_csv('${TRACE}', delim='\t', header=true, all_varchar=true, quote='');
CREATE TEMP TABLE agg AS
  SELECT name,
         count(*) FILTER (WHERE status = 'FAILED') AS attempts_failed,
         CASE WHEN count(*) FILTER (WHERE status IN ('COMPLETED','CACHED')) > 0
              THEN 'yes' ELSE 'no' END AS recovered
  FROM t GROUP BY name;
COPY (
  SELECT t.attempt, t.exit, t.name, t.tag, t.hash, t.workdir,
         regexp_replace(t.process, '^.*:', '') AS process,
         agg.attempts_failed, agg.recovered
  FROM t JOIN agg USING (name)
  WHERE t.status = 'FAILED'
  QUALIFY ROW_NUMBER() OVER (PARTITION BY t.name ORDER BY CAST(t.attempt AS INTEGER) DESC) = 1
  ORDER BY t.tag
) TO '${TMPD}/full.tsv' (DELIMITER '\t', HEADER true);
SQL
echo "==> duckdb: last failed attempt per task" >&2
duckdb -f "$TMPD/q.sql" >/dev/null
cut -f1-6 "$TMPD/full.tsv" > "$JOBS"
n_jobs=$(( $(wc -l < "$JOBS") - 1 ))
echo "    $n_jobs failed tasks -> $JOBS" >&2
if (( n_jobs == 0 )); then
    : > "$BLOB"
    printf 'process\ttag\texit\tattempts_failed\trecovered\tworkdir\tfirst_error\n' > "$MANI"
    echo "no failures in $RUN — nothing to triage" >&2
    exit 0
fi

# ---------------------------------------------------------------- 3. logs + manifest
# .command.err is the task's own stderr. .command.log under LSF has the whole
# LSF report and .command.run bolted on the end — that is what made the batch-1
# failed.log 227k lines for 1814 jobs and useless to grep. Truncate it.
first_error() {
    local f="$1" line=""
    # an explicit error line, wherever it is
    line="$(grep -m1 -E '\[ERROR\]|^ERROR\b|EXITING because of|Segmentation fault|Killed|TenxRunError|RenameError|^[[:space:]]*error:' "$f" 2>/dev/null || true)"
    # the exception line of a Python traceback is the LAST one
    [[ -z "$line" ]] && line="$(grep -E '^[A-Za-z_][A-Za-z0-9_.]*(Error|Exception|Exit):' "$f" 2>/dev/null | tail -1 || true)"
    # first real line, skipping the tools' own [INFO]/[WARN] progress chatter
    [[ -z "$line" ]] && line="$(grep -m1 -vE '^[[:space:]]*$|^\[(INFO|WARN|DEBUG)\]|^Traceback|^[[:space:]]+File "|^[[:space:]]+[~^]+[[:space:]]*$|^\+ ' "$f" 2>/dev/null || true)"
    # a job that only ever printed progress: its last line is where it died
    [[ -z "$line" ]] && line="$(grep -vE '^[[:space:]]*$' "$f" 2>/dev/null | tail -1 || true)"
    printf '%s' "${line//$'\t'/ }" | cut -c1-400
}

echo "==> building $BLOB and $MANI" >&2
: > "$BLOB"
printf 'process\ttag\texit\tattempts_failed\trecovered\tworkdir\tfirst_error\n' > "$MANI"
n_block=0; n_noblock=0; n_gone=0
while IFS=$'\t' read -r attempt exitcode name tag hash workdir process attempts_failed recovered; do
    src=""
    if   [[ -s "$workdir/.command.err" ]]; then src="$workdir/.command.err"
    elif [[ -s "$workdir/.command.log" ]]; then src="$workdir/.command.log"
    fi

    printf '### %s (%s) ###\n' "$tag" "$workdir" >> "$BLOB"
    if [[ -n "$src" ]]; then
        # stop at the LSF report; drop trailing blank / separator lines
        awk '/^Sender: LSF System/ {exit} {buf[NR]=$0}
             END {last=NR; while (last>0 && buf[last] ~ /^[-[:space:]]*$/) last--;
                  for (i=1; i<=last; i++) print buf[i]}' "$src" > "$TMPD/body"
        # exit 130 is an LSF kill, not a signal from the tool; the reason is
        # only ever in the report we just threw away.
        term="$(grep -m1 -oE 'TERM_[A-Z_]+' "$workdir/.command.log" 2>/dev/null || true)"
        [[ -n "$term" ]] && printf '[lsf] %s\n' "$term" >> "$TMPD/body"
        if [[ -s "$TMPD/body" ]]; then n_block=$((n_block+1)); else n_noblock=$((n_noblock+1)); fi
        cat "$TMPD/body" >> "$BLOB"
        fe="$(first_error "$TMPD/body")"
    else
        [[ -d "$workdir" ]] || n_gone=$((n_gone+1))
        n_noblock=$((n_noblock+1))
        echo "(no .command.log found)" >> "$BLOB"
        fe="(no .command.log found)"
    fi
    printf '\n' >> "$BLOB"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$process" "$tag" "$exitcode" "$attempts_failed" "$recovered" "$workdir" "$fe" >> "$MANI"
done < <(tail -n +2 "$TMPD/full.tsv")

# ---------------------------------------------------------------- 4. reconcile
n_mani=$(( $(wc -l < "$MANI") - 1 ))
n_recov=$(awk -F'\t' 'NR>1 && $5=="yes"' "$MANI" | wc -l)
cat >&2 <<EOF

run                : $RUN  (batch ${BATCH:-?}, suffix $SUFFIX)
failed tasks       : $n_jobs
  with stderr      : $n_block
  no log block     : $n_noblock$( (( n_gone )) && echo "   (work dir gone: $n_gone)" )
  recovered later  : $n_recov
manifest rows      : $n_mani   $( [[ $n_mani -eq $n_jobs ]] && echo OK || echo MISMATCH )

next: $(dirname "${BASH_SOURCE[0]}")/triage.py --manifest $MANI
EOF
[[ $n_mani -eq $n_jobs ]] || exit 1
