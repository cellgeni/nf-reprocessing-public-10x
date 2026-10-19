#!/bin/bash
# End-to-end check of nf-lineage-keep. Needs the plugin installed (make install)
# and Nextflow >= 26.04.1 on PATH. Runs everything locally in this directory.
#
#   1. a run in which ALIGN fails for good on S3
#   2. checks: S3's whole chain kept, every other chain emptied, scratch gone,
#      rerun.tsv = GSE2/S3
#   3. a -resume of rerun.tsv: S3's upstream tasks must come back CACHED
set -uo pipefail
cd "$(dirname "$0")"
rm -rf work results results-rerun .nextflow .nextflow.log*

printf 'dataset_id\tsample_id\nGSE1\t\nGSE2\tS3,S4\n' > datasets.tsv

fails=0
check() { if eval "$2"; then echo "ok   $1"; else echo "FAIL $1"; fails=$((fails+1)); fi; }

nextflow run main.nf -ansi-log false --datasets datasets.tsv --outdir results || { echo "run 1 failed"; exit 1; }

c=results/cleanup
echo; echo "== $c/rerun.tsv"; cat $c/rerun.tsv
echo "== $c/failed_tasks.tsv"; cat $c/failed_tasks.tsv
echo "== $c/deleted.tsv"; column -t -s $'\t' $c/deleted.tsv | cut -c1-200
echo "== $c/processes.tsv"; column -t -s $'\t' $c/processes.tsv

# work dir of a task, from the trace (last attempt)
wd() { awk -F'\t' -v n="$1" '$5==n {d=$9} END {print d}' results/trace.txt; }

check "rerun.tsv holds only GSE2/S3"        '[[ "$(tail -n +2 $c/rerun.tsv)" == $'"'"'GSE2\tS3'"'"' ]]'
check "failed_tasks.tsv lists ALIGN (S3)"   'grep -q "ALIGN (S3)" $c/failed_tasks.tsv'
for t in "DOWNLOAD (S3a)" "DOWNLOAD (S3b)" "RENAME_RUN (S3a)" "RENAME_RUN (S3b)" "RENAME_SAMPLE (S3)"; do
    check "kept: $t"                          '! grep -qF "$t" $c/deleted.tsv && [[ -n "$(ls $(wd "$t") | grep -v "^\." )" ]]'
done
for t in "DOWNLOAD (S1a)" "RENAME_RUN (S2b)" "RENAME_SAMPLE (S4)" "ALIGN (S1)" "QC (GSE1)" "QC (GSE2)"; do
    check "cleaned: $t"                       'grep -qF "$t" $c/deleted.tsv && [[ -z "$(ls -A $(wd "$t") | grep -v "^\.command\.\|^\.exitcode$")" ]]'
done
check "logs kept in a cleaned dir"          '[[ -f $(wd "ALIGN (S1)")/.command.sh ]]'
check "published outputs intact"            '[[ -s results/raw/GSE1/fastq/S1/S1a.fastq.gz && -s results/starsolo/GSE1/S1/matrix.mtx ]]'
check "scratch purged from failed ALIGN"    '! ls -d work/*/*/S3/_STARtmp 2>/dev/null | grep -q .'
check "failed ALIGN dir otherwise kept"     '[[ -f $(wd "ALIGN (S3)")/.command.err ]]'
# the downloaded FASTQ is hard-linked by publishing and by RENAME_RUN: counted, not freed
check "hard-linked bytes not counted freed" 'awk -F"\t" "\$2==\"delete\" && \$3==\"DOWNLOAD\" && (\$8 != 200000 || \$9 != 0) {bad=1} END {exit bad}" $c/deleted.tsv'

echo; echo "== resume with rerun.tsv"
nextflow run main.nf -ansi-log false --datasets $c/rerun.tsv --outdir results-rerun --fail none -resume || { echo "run 2 failed"; exit 1; }
status() { awk -F'\t' -v n="$1" '$5==n {s=$6} END {print s}' results-rerun/trace.txt; }
for t in "DOWNLOAD (S3a)" "DOWNLOAD (S3b)" "RENAME_RUN (S3a)" "RENAME_RUN (S3b)" "RENAME_SAMPLE (S3)"; do
    check "cached on resume: $t"              '[[ "$(status "$t")" == CACHED ]]'
done
check "ALIGN (S3) ran on resume"            '[[ "$(status "ALIGN (S3)")" == COMPLETED ]]'
check "nothing else ran on resume"          '! awk -F"\t" "NR>1" results-rerun/trace.txt | grep -qE "S1|S2|S4"'

echo; [[ $fails -eq 0 ]] && echo "ALL CHECKS PASSED" || { echo "$fails CHECK(S) FAILED"; exit 1; }
