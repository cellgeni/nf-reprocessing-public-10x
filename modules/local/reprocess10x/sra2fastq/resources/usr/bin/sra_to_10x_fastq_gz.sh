#!/bin/bash -e

mean_read_length() {
  local fq=$1
  head -n400 "$fq" | awk 'NR%4==2 {sum+=length($0); n++} END {if (n>0) printf "%d\n", sum/n+0.5; else print 0}'
}

merge_v1_cb_umi_fastq() {
  local cb_fq=$1
  local umi_fq=$2
  local out_fq=$3

  local cb_lines umi_lines
  cb_lines=$(wc -l < "$cb_fq")
  umi_lines=$(wc -l < "$umi_fq")

  if (( cb_lines % 4 != 0 || umi_lines % 4 != 0 )); then
    echo "ERROR: FASTQ line counts are invalid for v1 reconstruction: $cb_fq ($cb_lines), $umi_fq ($umi_lines)" >&2
    return 1
  fi

  if (( cb_lines != umi_lines )); then
    echo "ERROR: Barcode and UMI FASTQs have different numbers of records: $cb_fq ($cb_lines lines), $umi_fq ($umi_lines lines)" >&2
    return 1
  fi

  paste <(paste - - - - < "$cb_fq") <(paste - - - - < "$umi_fq") \
    | awk '
        BEGIN { FS=OFS="\t" }
        function norm_id(h,   x) {
          x=h
          sub(/^@/, "", x)
          sub(/[[:space:]].*$/, "", x)
          sub(/([\\/.][123])$/, "", x)
          return x
        }
        {
          cb_id  = norm_id($1)
          umi_id = norm_id($5)
          if (cb_id != umi_id) {
            printf("ERROR: barcode/UMI record ID mismatch at record %d: %s vs %s\n", NR, cb_id, umi_id) > "/dev/stderr"
            exit 1
          }
          print $1, $2 $6, $3, $4 $8
        }
      ' | tr '\t' '\n' > "$out_fq"
}

compress_all_fastqs() {
  local SRA=$1
  local CMD=${2:-""}
  local i
  for i in ${SRA}*fastq
  do
    [[ -e "$i" ]] || continue
    $CMD pigz "$i" &
  done
  wait
}

function sra2fastq {
  local SRA=$1
  local WL=$2
  local CPUS=$3
  local CMD=${4:-""}

  ## Identify 10x samples, reconstruct 3' v1 barcode+UMI reads for STARsolo,
  ## and otherwise keep standard barcode-in-read layouts unchanged.
  ##
  ## Important limitation:
  ##  - true Chromium 3' v1 can be detected from read geometry (14 bp CB in index + 10 bp UMI read)
  ##  - 5' GEX and 5' V(D)J both use barcode+UMI at the start of Read 1 and can share the same whitelist,
  ##    so they are not always distinguishable from FASTQ geometry alone
  ##  - therefore this script does not try to classify V(D)J; it only fixes 3' v1 and otherwise preserves
  ##    the standard barcode-in-read layout

  ## Dump FASTQ from SRA.
  ## NOTE: true 10x v1 requires technical/index reads to be present in the dump.
  $CMD parallel-fastq-dump -t $CPUS -T . --split-files -F -s $SRA

  ## Find which reads contain an actual 10x whitelist.
  local BC=""
  local MATCH_KIT=""
  local MATCH_COUNT=0
  local NMATCH=0
  local i N KIT

  for i in ${SRA}*fastq
  do
    [[ -e "$i" ]] || continue
    $CMD seqtk sample -s100 "$i" 200000 | awk 'NR%4==2' | cut -c-14 | grep -F -f "$WL/737K-april-2014_rc.txt" | wc -l > "$i.v1.count" &
    $CMD seqtk sample -s100 "$i" 200000 | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/737K-august-2016.txt"   | wc -l > "$i.v2.count" &
    $CMD seqtk sample -s100 "$i" 200000 | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/3M-february-2018.txt"   | wc -l > "$i.v3.count" &
    $CMD seqtk sample -s100 "$i" 200000 | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/737K-arc-v1.txt"        | wc -l > "$i.arc.count" &
    $CMD seqtk sample -s100 "$i" 200000 | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/3M-3pgex-may-2023.txt"  | wc -l > "$i.v4-3.count" &
    $CMD seqtk sample -s100 "$i" 200000 | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/3M-5pgex-jan-2023.txt"  | wc -l > "$i.v4-5.count" &
  done

  wait

  ## Analyze whitelist hits.
  for i in ${SRA}*count
  do
    [[ -e "$i" ]] || continue
    N=$(cat "$i")
    KIT=${i%%.count}
    KIT=${KIT##$SRA*.fastq.}
    if (( N > 50000 ))
    then
      BC=${i%%.$KIT.count}
      MATCH_KIT=$KIT
      MATCH_COUNT=$N
      echo "Barcode file: $BC, matching whitelist: $KIT, number of matched barcodes: $N"
      NMATCH=$((NMATCH+1))
    fi
  done

  ## Bail if more than 1 match; do not delete FASTQs if detection is ambiguous.
  if (( NMATCH > 1 ))
  then
    echo "WARNING: More than 1 file/whitelist match! This should not happen, please investigate the files in the log above."
  elif (( NMATCH == 0 ))
  then
    echo "WARNING: No files matched any of the whitelists! Most likely, this experiment is not a 10X single-cell RNA-seq."
  fi

  ## If one barcode read was identified, nominate the biological read.
  ## For 3' v1, also identify the 10 bp UMI read and reconstruct STARsolo-style barcode read = CB+UMI.
  if (( NMATCH == 1 ))
  then
    local BIO=""
    local BLEN=0
    local BCLEN=0
    local CURLEN=0
    local UMI=""
    local UMLEN=0
    local BESTDIFF=999
    local DIFF=0

    BCLEN=$(mean_read_length "$BC")

    for i in ${SRA}*fastq
    do
      [[ -e "$i" ]] || continue
      [[ "$i" == "$BC" ]] && continue
      CURLEN=$(mean_read_length "$i")
      if (( CURLEN >= BLEN ))
      then
        BIO=$i
        BLEN=$CURLEN
      fi
    done

    if [[ -z "$BIO" ]]
    then
      echo "WARNING: Could not identify a biological read after detecting barcode file $BC."
      compress_all_fastqs "$SRA" "$CMD"
      return 0
    fi

    echo "Detected barcode-like read: $BC (~${BCLEN} bp, whitelist=$MATCH_KIT, matches=$MATCH_COUNT)"
    echo "Longest remaining read chosen as biological candidate: $BIO (~${BLEN} bp)"

    ## Find a separate short UMI-like read if present.
    for i in ${SRA}*fastq
    do
      [[ -e "$i" ]] || continue
      [[ "$i" == "$BC" || "$i" == "$BIO" ]] && continue
      CURLEN=$(mean_read_length "$i")
      if (( CURLEN >= 9 && CURLEN <= 12 ))
      then
        DIFF=$(( CURLEN > 10 ? CURLEN - 10 : 10 - CURLEN ))
        if (( DIFF < BESTDIFF ))
        then
          UMI=$i
          UMLEN=$CURLEN
          BESTDIFF=$DIFF
        fi
      fi
    done

    if [[ "$MATCH_KIT" == "v1" ]]
    then
      if [[ -z "$UMI" ]]
      then
        echo "WARNING: Detected 10x 3' v1 barcode read ($BC) but could not find a 9-12 bp UMI FASTQ."
        echo "WARNING: This usually means the SRA dump did not preserve the technical/UMI read. Leaving all FASTQs compressed for manual inspection."
        compress_all_fastqs "$SRA" "$CMD"
        return 0
      fi

      echo "Detected official 10x Chromium 3' v1 layout: barcode=$BC (14 bp whitelist), UMI=$UMI (~${UMLEN} bp), cDNA=$BIO (~${BLEN} bp)"
      echo "Reconstructing STARsolo-compatible synthetic barcode read: 14 bp cell barcode + ${UMLEN} bp UMI ..."

      merge_v1_cb_umi_fastq "$BC" "$UMI" "${SRA}.tmp1"
      mv "$BIO" "${SRA}.tmp2"

      if ls ${SRA}*fastq >/dev/null 2>&1
      then
        rm -f ${SRA}*fastq
      fi

      mv "${SRA}.tmp1" "${SRA}_1.fastq"
      mv "${SRA}.tmp2" "${SRA}_2.fastq"

      $CMD pigz "${SRA}_1.fastq" &
      $CMD pigz "${SRA}_2.fastq" &
    else
      ## Standard barcode-in-read layouts.
      ## This covers 3' v2/v3/v4 and also 5' chemistries, whose raw FASTQs already provide CB+UMI in one read.
      ## The script intentionally does not try to distinguish 5' GEX from 5' V(D)J here.
      if (( BCLEN >= 27 && BCLEN <= 30 && BLEN >= 70 )); then
        echo "Detected standard barcode-in-read layout consistent with 5' or late 3' chemistries; keeping reads as-is."
      else
        echo "Detected standard barcode-in-read layout; keeping reads as-is."
        echo "NOTE: 737K-august-2016 is shared by 3' v2 and 5' v1/v2, so raw FASTQ geometry alone cannot always separate them."
      fi

      mv "$BC" "${SRA}.tmp1"
      mv "$BIO" "${SRA}.tmp2"
      if ls ${SRA}*fastq >/dev/null 2>&1
      then
        rm -f ${SRA}*fastq
      fi
      mv "${SRA}.tmp1" "${SRA}_1.fastq"
      mv "${SRA}.tmp2" "${SRA}_2.fastq"

      $CMD pigz "${SRA}_1.fastq" &
      $CMD pigz "${SRA}_2.fastq" &
    fi
  else
    ## If things look weird, just compress everything, we'll sort it later.
    compress_all_fastqs "$SRA" "$CMD"
  fi

  wait
}

function main () {
  PARSED=$1
  SRA=$(grep -w "SRA$" "$PARSED" | cut -f1 | head -"$LSB_JOBINDEX" | tail -1)
  SIF="/nfs/cellgeni/singularity/images/reprocess_10x.sif"
  CMD="singularity run --bind /nfs,/lustre $SIF"
  WL=/nfs/cellgeni/STAR/whitelists
  CPUS=16

  sra2fastq "$SRA" "$WL" "$CPUS" "$CMD"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
