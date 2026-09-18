#!/bin/bash
# Post a "pipeline finished" ping to Slack, with what is left of the Lustre
# group quota. Meant to be the line right after `nextflow run` in
# scripts/run_reprocess.bsub, so a finished (or dead) run announces itself
# instead of being noticed hours later:
#
#   nextflow run main.nf ... ; status=$?
#   .claude/skills/reprocess-debug/bin/notify_run_done.sh --batch "$batch" --exit-code "$status"
#   exit $status
#
# This is NOT the debug notification. notify.py + references/notify.md are for
# the end of a debug session (post-mortem link, failure attribution, project
# progress); this is a bare "it stopped, here is how much scratch is left",
# with no triage in it, and it must stay cheap and never fail the job: every
# lookup here degrades to "unknown" rather than exiting non-zero.
#
# Quota is read with `lfs quota -q -g <group> <fs>`, whose numbers are KiB and
# whose columns are  used soft hard grace files soft hard grace  after the
# filesystem path. The soft limit is what actually bites first, so it is
# preferred over the hard limit when it is set (on scratch124 it is 0 = unset,
# and the 142T hard limit is the real ceiling). A `*` suffix on a used value
# means "over soft quota, in grace" — it is stripped before arithmetic.
set -uo pipefail

FS="/lustre/scratch124/cellgen"
GROUP="cellgeni"
SLACK_ENV_FILE="/lustre/scratch124/cellgen/cellgeni/aljes/reprocessing_slack_bot/.env"
BATCH="" ; RUN="" ; EXIT_CODE="" ; TO="" ; DRY_RUN=0
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat >&2 <<'EOF'
usage: notify_run_done.sh [options]
  --batch <name|N>      batch label for the message ("batch10" or "10")
  --run <name>          Nextflow run name; default: last row of .nextflow/history
  --exit-code <n>       nextflow's exit status ($? of the run)
  --to <addr>           also send the same thing as an email (default: Slack only)
  --fs <path>           filesystem to report quota for (default: /lustre/scratch124/cellgen)
  --group <name>        quota group (default: cellgeni)
  --slack-env-file <f>  KEY=VALUE file with SLACK_BOT_TOKEN/CHANNEL_ID
  --dry-run             print the message, send nothing
EOF
}

while (( $# )); do
    case "$1" in
        --batch)          BATCH="$2"; shift 2 ;;
        --run)            RUN="$2"; shift 2 ;;
        --exit-code)      EXIT_CODE="$2"; shift 2 ;;
        --to)             TO="$2"; shift 2 ;;
        --fs)             FS="$2"; shift 2 ;;
        --group)          GROUP="$2"; shift 2 ;;
        --slack-env-file) SLACK_ENV_FILE="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "notify_run_done: unknown option '$1'" >&2; usage; exit 2 ;;
    esac
done

[[ -n "$BATCH" && "$BATCH" =~ ^[0-9]+$ ]] && BATCH="batch${BATCH}"

# --- quota ------------------------------------------------------------------
# One line: used_h limit_h free_h used_pct files_used_h files_limit_h files_pct
quota_fields() {
    command -v lfs >/dev/null 2>&1 || return 1
    lfs quota -q -g "$GROUP" "$FS" 2>/dev/null \
      | tr -s '[:space:]' '\n' | grep -v '^$' \
      | awk '
        function h(k,   u, i) {          # KiB -> human, same units as lfs -h
            split("K M G T P", u, " "); i = 1
            while (k >= 1024 && i < 5) { k /= 1024; i++ }
            return sprintf("%.1f%s", k, u[i])
        }
        function hn(n,   u, i) {          # bare count -> human
            # "|" and not " ": awk`s split() on a single space strips the empty
            # leading field, which silently shifts every unit up by one.
            split("|k|M|G", u, "|"); i = 1
            while (n >= 1000 && i < 4) { n /= 1000; i++ }
            return (i == 1) ? sprintf("%d", n) : sprintf("%.1f%s", n, u[i])
        }
        $1 ~ /^\// { seen = 1; n = 0; next }
        seen && n < 8 { t = $1; sub(/\*$/, "", t); v[++n] = t }
        END {
            if (n < 8) exit 1
            used = v[1] + 0; cap = (v[2] + 0 > 0) ? v[2] + 0 : v[3] + 0
            f_used = v[5] + 0; f_cap = (v[6] + 0 > 0) ? v[6] + 0 : v[7] + 0
            printf "%s %s %s %s %s %s %s\n",
                h(used),
                (cap > 0 ? h(cap) : "unlimited"),
                (cap > 0 ? h(cap - used) : "unlimited"),
                (cap > 0 ? sprintf("%.0f%%", 100 * used / cap) : "n/a"),
                hn(f_used),
                (f_cap > 0 ? hn(f_cap) : "unlimited"),
                (f_cap > 0 ? sprintf("%.0f%%", 100 * f_used / f_cap) : "n/a")
        }'
}

Q_USED="?" ; Q_CAP="?" ; Q_FREE="?" ; Q_PCT="?" ; Q_FILES="?" ; Q_FCAP="?" ; Q_FPCT="?"
if read -r a b c d e f g < <(quota_fields) && [[ -n "${g:-}" ]]; then
    Q_USED="$a"; Q_CAP="$b"; Q_FREE="$c"; Q_PCT="$d"; Q_FILES="$e"; Q_FCAP="$f"; Q_FPCT="$g"
else
    echo "notify_run_done: could not read 'lfs quota -g $GROUP $FS' — reporting unknown" >&2
fi

# --- run name, duration, nextflow's own verdict -----------------------------
# .nextflow/history is timestamp \t duration \t run_name \t status \t ... ; the
# duration and status columns stay "-" for a run that was killed rather than
# finishing, which is exactly the case worth being told about.
STARTED="unknown" ; DURATION="unknown" ; NF_STATUS="unknown"
if [[ -f .nextflow/history ]]; then
    row=""
    if [[ -n "$RUN" ]]; then
        row="$(awk -F'\t' -v r="$RUN" '$3 == r {line=$0} END {print line}' .nextflow/history)"
    elif [[ -n "$BATCH" ]]; then
        row="$(awk -F'\t' -v b="${BATCH}.csv" 'index($7, b) {line=$0} END {print line}' .nextflow/history)"
    fi
    [[ -z "$row" ]] && row="$(tail -1 .nextflow/history)"
    if [[ -n "$row" ]]; then
        STARTED="$(cut -f1 <<<"$row")"
        DURATION="$(cut -f2 <<<"$row")"
        RUN="${RUN:-$(cut -f3 <<<"$row")}"
        NF_STATUS="$(cut -f4 <<<"$row")"
        [[ "$DURATION"  == "-" ]] && DURATION="unknown (no end recorded — killed?)"
        [[ "$NF_STATUS" == "-" ]] && NF_STATUS="no status recorded (killed?)"
    fi
fi
RUN="${RUN:-unknown}"
BATCH="${BATCH:-unknown batch}"

if [[ -z "$EXIT_CODE" ]]; then
    VERDICT_EMOJI="🏁" ; VERDICT="finished — exit status not passed in"
elif [[ "$EXIT_CODE" == "0" ]]; then
    VERDICT_EMOJI="✅" ; VERDICT="finished, exit 0"
else
    VERDICT_EMOJI="❌" ; VERDICT="exited $EXIT_CODE"
fi

# errorStrategy 'ignore' keeps task failures out of the exit status (CLAUDE.md),
# so exit 0 says nothing about how many samples died. Say that here rather than
# letting a green tick read as "nothing failed".
CAVEAT="exit status ignores per-task failures (errorStrategy 'ignore') — triage the run before trusting it"

SLACK_FILE="$(mktemp -t notify-run-done-slack.XXXXXX)"
MAIL_FILE="$(mktemp -t notify-run-done-mail.XXXXXX)"
trap 'rm -f "$SLACK_FILE" "$MAIL_FILE"' EXIT

cat >"$SLACK_FILE" <<EOF
${VERDICT_EMOJI} *Pipeline run ${VERDICT}* — \`${BATCH}\` / \`${RUN}\`
🕐 Started ${STARTED} · ⏱ ${DURATION} · 📄 Nextflow says: ${NF_STATUS}
🖥 \`${HOSTNAME:-$(hostname)}\`${LSB_JOBID:+ · LSF job \`$LSB_JOBID\`}

💾 *Lustre* \`${FS}\` (group \`${GROUP}\`)
📉 Used: *${Q_USED}* / ${Q_CAP} (${Q_PCT}) — *${Q_FREE} free*
🗂 Files: ${Q_FILES} / ${Q_FCAP} (${Q_FPCT})

⚠️ _${CAVEAT}_
EOF

cat >"$MAIL_FILE" <<EOF
${VERDICT_EMOJI} PIPELINE RUN ${VERDICT^^}
================================================================

Batch:     ${BATCH}
Run:       ${RUN}
Started:   ${STARTED}
Duration:  ${DURATION}
Nextflow:  ${NF_STATUS}
Host:      ${HOSTNAME:-$(hostname)}${LSB_JOBID:+
LSF job:   $LSB_JOBID}

💾 LUSTRE ${FS} (group ${GROUP})
- Used:      ${Q_USED} / ${Q_CAP} (${Q_PCT})
- Free:      ${Q_FREE}
- Files:     ${Q_FILES} / ${Q_FCAP} (${Q_FPCT})

⚠️ Caveat: ${CAVEAT}

--
sent by notify_run_done.sh from scripts/run_reprocess.bsub 🤖
EOF

if (( DRY_RUN )); then
    echo "--- Slack ---" ; cat "$SLACK_FILE"
    [[ -n "$TO" ]] && { echo "--- email to $TO ---" ; cat "$MAIL_FILE" ; }
    exit 0
fi

# notify.py owns the credential handling and the uvx fallback for slack_sdk; a
# second copy of either would drift. It exits 0 when Slack is skipped, so a
# missing token never fails the LSF job either.
args=(--slack --slack-text-file "$SLACK_FILE" --slack-env-file "$SLACK_ENV_FILE")
if [[ -n "$TO" ]]; then
    args+=(--to "$TO" --subject "${BATCH}: pipeline run ${VERDICT} (${Q_FREE} free on scratch124)"
           --body-file "$MAIL_FILE")
fi
"$HERE/notify.py" "${args[@]}" || echo "notify_run_done: notify.py failed — not failing the job over it" >&2
exit 0
