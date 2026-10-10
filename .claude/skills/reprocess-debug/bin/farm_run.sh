#!/bin/bash
# Run one measurement on an LSF execution node instead of the head node.
#
#   farm_run.sh --name rle-SRR17720155 -- awk 'NR%4==2{...}' nf-work/65/c3fc80…/SRR17720155_2.fastq
#   farm_run.sh --name matstats --mem 8G -- python3 /path/to/matstats.py cands.txt out.json
#   farm_run.sh --name collect23 -- .claude/skills/reprocess-debug/bin/collect_run_logs.sh --batch 23
#   farm_run.sh --no-wait --name gzt -- 'for f in nf-work/3e/…/*.gz; do gzip -t "$f" && echo "$f OK"; done'
#
# The head nodes are for editing, LSF and light reads. On 2026-10-07 a batch-23
# debug session ran full-file awk scans of 30-110 GB FASTQs, gzip -t and a
# python pass over 218 matrices there, in the background, and Arbiter put ab76
# in penalty1 (CPU cut to 80% for 30 min; awk alone averaged ~300% CPU).
# "In the background" is not "off the head node".
#
# Submits with bsub -K, which blocks until the job ends and returns its exit
# status, so from an agent session run this with run_in_background: true and
# you are told when it finishes. The command's own stdout is printed at the end
# (tail only if it is large); stdout, stderr, the LSF report and bsub's own
# chatter are kept in --outdir (default logs/lsf/debug/) as
# <name>.<stamp>.{out,err,lsf,bsub}.
#
# One argument after -- is run as a bash snippet (pipes, loops, globs);
# several are run as an argv. The job starts in the current directory and
# inherits this shell's environment, modules included. It cannot see the head
# node's /tmp (the agent scratchpad lives there): keep helper scripts and
# outputs on Lustre, e.g. in logs/lsf/debug/. The wrapper refuses /tmp paths.
set -uo pipefail

NAME="job"; MEM="4G"; CPUS=1; QUEUE="normal"; WALL="4:00"; GROUP="cellgeni"
WAIT=1; DRY=0; ALLOW_TMP=0
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
OUTDIR="$REPO/logs/lsf/debug"
PRINT_MAX=$((256 * 1024))

die() { echo "error: $*" >&2; exit 1; }

usage() {
    awk 'NR==1{next} /^#/{sub(/^# ?/,""); print; next} {exit}' "${BASH_SOURCE[0]}"
    echo
    echo "options: [--name N] [--mem 4G] [--cpus 1] [--queue normal] [--time H:MM]"
    echo "         [--group cellgeni] [--outdir DIR] [--no-wait] [--dry-run] [--allow-tmp] -- command…"
    echo "queues : normal (12 h, default), long (48 h), transfer (outbound network)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)    NAME="$2";   shift 2 ;;
        --mem)     MEM="$2";    shift 2 ;;
        --cpus)    CPUS="$2";   shift 2 ;;
        --queue)   QUEUE="$2";  shift 2 ;;
        --time)    WALL="$2";   shift 2 ;;
        --group)   GROUP="$2";  shift 2 ;;
        --outdir)  OUTDIR="$2"; shift 2 ;;
        --no-wait) WAIT=0; shift ;;
        --dry-run) DRY=1;  shift ;;
        --allow-tmp) ALLOW_TMP=1; shift ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) die "unknown argument: $1 (commands go after --; try --help)" ;;
    esac
done
[[ $# -gt 0 ]] || usage
command -v bsub >/dev/null 2>&1 || die "bsub not on PATH — this needs an LSF submit host"

# LSF wants -M and rusage in MB here; accept 4G / 4000M / 4000
case "$MEM" in
    *[Gg]) MEM_MB=$(( ${MEM%[Gg]} * 1024 )) ;;
    *[Mm]) MEM_MB=${MEM%[Mm]} ;;
    *[0-9]) MEM_MB=$MEM ;;
    *) die "--mem: expected e.g. 4G or 4000M, got '$MEM'" ;;
esac
[[ $NAME =~ ^[A-Za-z0-9._-]+$ ]] || die "--name: letters, digits, . _ - only"

# /tmp is node-local. A job reading a helper script from the session scratchpad
# finds nothing, and one writing there writes to the execution node's /tmp,
# where the result is lost (batch 23's first LSF collection did exactly that).
# Keep job inputs and outputs on Lustre/NFS, e.g. logs/lsf/debug/.
case "$PWD/" in /tmp/*) die "the job would start in node-local $PWD; cd to the repo or Lustre first" ;; esac
case "$OUTDIR/" in /tmp/*) die "--outdir $OUTDIR is node-local; use a Lustre/NFS path" ;; esac
for a in "$@"; do
    [[ $a == *"/tmp/"* ]] && (( ! ALLOW_TMP )) && die "the command mentions a /tmp path, which an execution node cannot see:
       $a
       put helper scripts and outputs under logs/lsf/debug/ (or pass --allow-tmp if the job only uses its own /tmp)"
done

mkdir -p "$OUTDIR" || die "cannot create $OUTDIR"
STAMP=$(date +%Y%m%d-%H%M%S)
SCRIPT="$OUTDIR/$NAME.$STAMP.sh"

# The job script writes the command's streams to files of its own, so the
# LSF report (which LSF appends to -o) never mixes into the output.
{
    echo '#!/bin/bash'
    echo "cd $(printf '%q' "$PWD")"
    echo 'out="$1"; err="$2"'
    if [[ $# -eq 1 ]]; then
        printf 'bash -c %q > "$out" 2> "$err"\n' "$1"
    else
        printf '%q ' "$@"; echo '> "$out" 2> "$err"'
    fi
} > "$SCRIPT"
chmod +x "$SCRIPT"

PFX="$OUTDIR/$NAME.$STAMP"
BSUB=(bsub -G "$GROUP" -q "$QUEUE" -n "$CPUS" -W "$WALL"
      -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB] span[hosts=1]"
      -J "rdebug-$NAME" -o "$PFX.lsf" -e "$PFX.lsf")
(( WAIT )) && BSUB+=(-K)
BSUB+=("$SCRIPT" "$PFX.out" "$PFX.err")

if (( DRY )); then
    printf '%q ' "${BSUB[@]}"; echo
    echo "--- $SCRIPT"; cat "$SCRIPT"
    rm -f "$SCRIPT"
    exit 0
fi

echo "submitting rdebug-$NAME: queue $QUEUE, $CPUS cpu, $MEM, -W $WALL" >&2
if (( ! WAIT )); then
    "${BSUB[@]}"
    status=$?
    echo "outputs: $PFX.{out,err,lsf}   (bjobs -J rdebug-$NAME to follow)" >&2
    exit $status
fi

"${BSUB[@]}" > "$PFX.bsub" 2>&1
status=$?
grep -oE 'Job <[0-9]+>' "$PFX.bsub" | head -1 >&2

reason=$(grep -oE 'TERM_[A-Z_]+' "$PFX.lsf" 2>/dev/null | head -1)
[[ -n $reason ]] && echo "lsf: $reason" >&2
if [[ -s $PFX.err ]]; then
    echo "--- stderr (last 20 lines, $PFX.err)" >&2
    tail -n 20 "$PFX.err" >&2
fi
if [[ -f $PFX.out ]]; then
    size=$(stat -c %s "$PFX.out")
    if (( size > PRINT_MAX )); then
        echo "--- stdout is $size bytes; last 200 lines (full: $PFX.out)" >&2
        tail -n 200 "$PFX.out"
    else
        cat "$PFX.out"
    fi
else
    echo "no stdout file — the job never started? see $PFX.bsub and $PFX.lsf" >&2
fi
echo "exit $status · outputs $PFX.{out,err,lsf}" >&2
exit $status
