#!/usr/bin/env bash
# Regression test for GitHub #41: two BAMs of one sample covering the same lane must not
# produce FASTQs with the same name.
#
# REPROCESS10X_BAM2FASTQ renames each BAM's bamtofastq output with rename_fastqs.sh, and
# RENAME10XSAMPLE then stages every FASTQ of the sample into one flat directory
# (`stageAs: "fastqs/*"`). Two identical basenames there are an input file name collision,
# which Nextflow raises before the task runs, so it aborts the whole run instead of failing
# one sample. This test rebuilds that path without Nextflow or data:
#
#   1. fake bamtofastq layouts for three runs of one sample, all on lane L001 (one run with
#      two bamtofastq subdirectories, to cover rename_fastqs.sh's S-counter);
#   2. rename_fastqs.sh on each, with the prefix argument main.nf actually passes;
#   3. every output linked into one flat directory, as RENAME10XSAMPLE stages them;
#   4. rename_fastqs_recommended.py on that directory, which must give one lane per run.
#
# Usage: bash modules/local/reprocess10x/bam2fastq/tests/rename_fastqs.test.sh
#        (from anywhere; needs only bash and python3, so it is light enough for a head node)

set -euo pipefail

MODULE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
REPO=$(cd "$MODULE/../../../.." && pwd)
RENAME_RUN="$MODULE/resources/usr/bin/rename_fastqs.sh"
RENAME_SAMPLE="$REPO/modules/cellgeni/rename10xsample/resources/usr/bin/rename_fastqs_recommended.py"
SAMPLE=ERS3741160

fail() { echo "FAIL: $*" >&2; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Which meta key main.nf passes as rename_fastqs.sh's prefix. Reading it from the module,
# rather than restating it here, is what makes this test fail if the call is reverted.
prefix_key=$(sed -n 's/^[[:space:]]*rename_fastqs\.sh \${meta\.\([a-z_]*\)} \${meta\.id} fastqs.*/\1/p' "$MODULE/main.nf")
case "$prefix_key" in
    id|sample_id) ;;
    *) fail "could not find 'rename_fastqs.sh \${meta.<key>} \${meta.id} fastqs' in $MODULE/main.nf" ;;
esac

# A small but valid gzipped FASTQ: rename_fastqs_recommended.py samples read lengths from it.
make_fastq() {  # <path> <read length>
    local seq qual i
    seq=$(printf 'A%.0s' $(seq 1 "$2"))
    qual=$(printf 'F%.0s' $(seq 1 "$2"))
    for i in 1 2 3 4; do printf '@read%d\n%s\n+\n%s\n' "$i" "$seq" "$qual"; done | gzip > "$1"
}

# bamtofastq writes <run>/<run>_<gem>_<lib>_<flowcell>/bamtofastq_S1_L00n_<read>_001.fastq.gz,
# with S1 in every subdirectory. rename_fastqs.sh runs from the task dir with <run> as source.
make_run() {  # <run> <number of bamtofastq subdirectories>
    local run=$1 n=$2 k dir
    for k in $(seq 1 "$n"); do
        dir="$WORK/$run/$run/${run}_0_${k}_HXXXXXXX${k}"
        mkdir -p "$dir"
        make_fastq "$dir/bamtofastq_S1_L001_I1_001.fastq.gz" 8
        make_fastq "$dir/bamtofastq_S1_L001_R1_001.fastq.gz" 28
        make_fastq "$dir/bamtofastq_S1_L001_R2_001.fastq.gz" 90
    done
    local prefix
    if [[ $prefix_key == id ]]; then prefix=$run; else prefix=$SAMPLE; fi
    (cd "$WORK/$run" && bash "$RENAME_RUN" "$prefix" "$run" fastqs > rename.log)
}

# ERR3517557 and ERR3517561 are SA220_S4_L001.bam and SA220_S7_L001.bam, the pair behind the
# REQ-81750 abort. ERR0000001 is a made-up BAM holding two libraries.
make_run ERR3517557 1
make_run ERR3517561 1
make_run ERR0000001 2

# 3. Stage like RENAME10XSAMPLE: one flat directory. A duplicate name is the #41 collision.
dups=$(find "$WORK"/ERR*/fastqs -name '*.fastq.gz' -printf '%f\n' | sort | uniq -d)
[[ -z $dups ]] || fail "BAM2FASTQ outputs collide when staged together (prefix: meta.$prefix_key):"$'\n'"$dups"
mkdir "$WORK/staged"
ln "$WORK"/ERR*/fastqs/*.fastq.gz "$WORK/staged/"

# 4. The sample-level rename must accept them and give each run (and each subdirectory) a lane.
mkdir "$WORK/out"
python3 "$RENAME_SAMPLE" --fastqs "$WORK"/staged/*.fastq.gz --sample-id "$SAMPLE" \
    --outdir "$WORK/out" 2> "$WORK/rename_sample.log" \
    || { cat "$WORK/rename_sample.log" >&2; fail "rename_fastqs_recommended.py exited non-zero"; }

expected=$(for lane in 1 2 3 4; do for read in I1 R1 R2; do
    printf '%s_S1_L%03d_%s_001.fastq.gz\n' "$SAMPLE" "$lane" "$read"; done; done)
actual=$(ls "$WORK/out")
[[ $actual == "$expected" ]] || fail "unexpected sample-level outputs:"$'\n'"$actual"$'\n'"expected:"$'\n'"$expected"

echo "PASS: 3 BAMs (one with two libraries), all on L001, staged without collision and renamed to L001-L004"
