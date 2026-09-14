#!/bin/bash
# Sourced helper. Resolving a batch number to a Nextflow run name is shared by
# collect_run_logs.sh and archive_run.sh, and it is subtle enough that two
# copies would drift: a -resume keeps the session uuid and takes a NEW run name,
# so a batch routinely has more than one run, and picking the wrong one either
# shows failures a later resume fixed, or hides them.
#
#   resolve_run <batch> <run> <latest>
#     -> RR_RUN RR_BATCH RR_ROW   (RR_ROW is the raw .nextflow/history line)
#
# Exactly one of <batch>/<run> need be set. <latest>=1 takes the most recent run
# of a batch instead of refusing to choose. Call from the repo root.

# .nextflow/history: timestamp \t duration \t run_name \t status \t script_id \t session \t command
resolve_run() {
    local batch="$1" run="$2" latest="${3:-0}" row="" ; local -a cands

    [[ -f .nextflow/history ]] || { echo "error: no .nextflow/history here" >&2; return 1; }

    if [[ -z "$run" ]]; then
        [[ -n "$batch" ]] || { echo "error: resolve_run needs a batch or a run" >&2; return 1; }
        mapfile -t cands < <(awk -F'\t' -v b="batch${batch}.csv" 'index($7, b) {print}' .nextflow/history)
        (( ${#cands[@]} )) || { echo "error: no run in .nextflow/history mentions batch${batch}.csv" >&2; return 1; }
        if (( ${#cands[@]} > 1 && latest == 0 )); then
            echo "batch $batch has ${#cands[@]} runs (a -resume keeps the session and takes a new name):" >&2
            printf '%s\n' "${cands[@]}" | awk -F'\t' '{printf "  %s  %-22s %-4s %s\n", $1, $3, $4, $2}' >&2
            echo "error: pick one with --run <name>, or take the most recent with --latest" >&2
            return 1
        fi
        row="${cands[-1]}"
        run="$(cut -f3 <<<"$row")"
    else
        row="$(awk -F'\t' -v r="$run" '$3 == r {print}' .nextflow/history | tail -1)"
        [[ -n "$row" ]] || { echo "error: run '$run' not in .nextflow/history" >&2; return 1; }
        [[ -z "$batch" ]] && batch="$(grep -oE 'batch[0-9]+\.csv' <<<"$row" | head -1 | grep -oE '[0-9]+' || true)"
    fi

    RR_RUN="$run"; RR_BATCH="$batch"; RR_ROW="$row"
}

# The LSF driver job id is a third identity for the same run: logs/lsf/
# reprocessOutput<JOBID>.log carries "Launching `main.nf` [<run_name>]" near the top.
resolve_lsf_jobid() {
    local run="$1" f id
    for f in logs/lsf/reprocessOutput*.log; do
        [[ -e "$f" ]] || continue
        if grep -q "\[${run}\]" "$f" 2>/dev/null; then
            id="$(basename "$f" | tr -dc '0-9')"; echo "$id"; return 0
        fi
    done
    return 1
}

# Nextflow names its reports after the run's start time, but NOT after the
# timestamp .nextflow/history records — they differ by a second or three, and
# once by 3s (peaceful_feynman: history 13:27:57, report 13:27:54). Matching the
# stamp exactly finds a report for one run in five and silently reports the rest
# as missing. So match the nearest stamp within a tolerance instead.
#
#   find_report_stamp "<history timestamp>" [tolerance_seconds]
#     -> prints the stamp of the closest reports/execution_report_*.html, or
#        nothing (exit 1) if none is within tolerance.
find_report_stamp() {
    local ts="${1%%$'\t'*}" tol="${2:-300}" want f st best="" bestd=""
    want=$(date -d "$ts" +%s 2>/dev/null) || return 1
    for f in reports/execution_report_*.html; do
        [[ -e "$f" ]] || continue
        st="${f##*/execution_report_}"; st="${st%.html}"
        # "2026-09-08_12-43-36" -> "2026-09-08 12:43:36"
        local human="${st%%_*} ${st#*_}"; human="${human//-/:}"
        human="${human:0:10} ${human:11}"
        human="$(echo "$st" | sed 's/_/ /; s/\(.*\) \([0-9]*\)-\([0-9]*\)-\([0-9]*\)/\1 \2:\3:\4/')"
        local e d
        e=$(date -d "$human" +%s 2>/dev/null) || continue
        d=$(( want > e ? want - e : e - want ))
        if [[ -z "$bestd" || $d -lt $bestd ]]; then bestd=$d; best="$st"; fi
    done
    [[ -n "$best" && $bestd -le $tol ]] || return 1
    printf '%s' "$best"
}
