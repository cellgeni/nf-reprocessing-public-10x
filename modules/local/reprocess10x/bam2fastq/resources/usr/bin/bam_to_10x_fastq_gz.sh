#!/bin/bash


## Cell Ranger >= 1.2 and Long Ranger >= 2.1 stamp '@CO 10x_bam_to_fastq:' lines into the
## BAM header telling bamtofastq which tag holds which read. Some archived BAMs — notably
## the pre-alignment ones written by 'fastq_pre_barcodes', which several ArrayExpress
## submissions publish instead of FASTQs — carry every piece of the read in tags but
## describe none of it. bamtofastq then reports "Unrecognized 10x BAM file" and suggests
## --gemcode / --lr20 / --cr11. All three are wrong for these: those flags assume a 14 bp
## v1 barcode, while these BAMs hold a 16 bp barcode in CR with the UMI in RX.
##
## Everything bamtofastq needs is present, just undeclared:
##   CR:Z / CY:Z   cell barcode + quality      -> R1, together with the UMI
##   RX:Z / QX:Z   UMI + quality               (UR:UY in Cell Ranger's own naming)
##   BC:Z / QT:Z   sample index + quality      -> I1
##   SEQ  / QUAL   biological read             -> R2
##
## Two conditions both have to hold before bamtofastq emits anything, and each fails
## silently on its own:
##   1. the header carries '@CO 10x_bam_to_fastq:' lines describing the layout; without
##      them the file is "unrecognized" and nothing is written;
##   2. every record carries an RG:Z tag naming an @RG whose ID is Cell Ranger's
##      five-field sample:library:gem_group:flowcell:lane form. An @RG with a plain ID is
##      ignored outright ("no @RG headers found"), and with no usable read group
##      bamtofastq falls back to splitting on the corrected-barcode tag CB — which these
##      BAMs do not have. It then reads every record and writes zero.
##
## So the patch below both reheaders the file and stamps RG:Z onto every record. That is a
## full rewrite of the BAM; it is done in a single streaming pass and the copy is deleted
## as soon as the FASTQs exist, keeping peak disk at roughly twice the input.

## Does the header already describe its own layout? If so the BAM needs no help.
function header_describes_10x_layout {
  local HEADER=$1
  grep -q '^@CO[[:space:]]*10x_bam_to_fastq:' <<< "$HEADER"
}

## Rewrite $RUN.bam as $RUN.patched.bam with a layout the vendor tool understands.
## Fails loudly if the tags a 10x GEX read needs are not all there — an honest error
## beats bamtofastq's misleading "use --cr11" suggestion.
function patch_legacy_10x_bam {
  local RUN=$1
  local CPUS=$2
  local CMD=${3:-""}
  local HEADER=$4

  ## One record is enough to read the tag layout; SIGPIPE from head is expected here.
  local FIRST
  FIRST=$($CMD samtools view "$RUN.bam" 2>/dev/null | head -1) || true

  if [[ -z $FIRST ]]
  then
    >&2 echo "ERROR: $RUN.bam has no 10x_bam_to_fastq header comments and no records to inspect."
    return 1
  fi

  if ! grep -q $'\tCR:Z:' <<< "$FIRST" || ! grep -q $'\tCY:Z:' <<< "$FIRST"
  then
    >&2 echo "ERROR: $RUN.bam is not a 10x BAM this pipeline can convert: it has neither"
    >&2 echo "       '@CO 10x_bam_to_fastq:' header comments nor a CR/CY cell-barcode tag."
    return 1
  fi

  ## Cell Ranger writes the UMI as UR/UY; fastq_pre_barcodes writes it as RX/QX.
  local UMI_TAG UMI_QUAL
  if grep -q $'\tRX:Z:' <<< "$FIRST"
  then
    UMI_TAG=RX; UMI_QUAL=QX
  elif grep -q $'\tUR:Z:' <<< "$FIRST"
  then
    UMI_TAG=UR; UMI_QUAL=UY
  else
    >&2 echo "ERROR: $RUN.bam has a CR cell-barcode tag but no UMI tag (looked for RX and UR)."
    return 1
  fi

  ## The sample index is optional — plenty of submissions drop I1 — so only declare it
  ## when it is actually there. Declaring it otherwise makes bamtofastq drop every read.
  local HAVE_I1=no
  if grep -q $'\tBC:Z:' <<< "$FIRST" && grep -q $'\tQT:Z:' <<< "$FIRST"
  then
    HAVE_I1=yes
  fi

  ## bamtofastq takes the lane from the 5th field of the @RG ID and puts it in the output
  ## filename, so it is worth getting right. The @PG line of a fastq_pre_barcodes BAM
  ## names the FASTQs it was built from (..._S1_L002_R2_001.fastq.gz), which is the most
  ## reliable lane we have; anything else falls back to lane 1. The flowcell field only
  ## names an output subdirectory that rename_fastqs.sh discards, so the run id serves.
  local LANE
  LANE=$(grep -o '_L[0-9][0-9][0-9]_' <<< "$HEADER" | head -1 | tr -dc '0-9' | sed 's/^0*//') || true
  if [[ -z ${LANE:-} ]]
  then
    LANE=1
  fi

  local RGID="$RUN:0:1:$RUN:$LANE"

  >&2 echo "NOTE: $RUN.bam carries no 10x_bam_to_fastq header comments — treating it as a"
  >&2 echo "      legacy pre-barcode BAM (CB=CR:CY, UMI=$UMI_TAG:$UMI_QUAL, I1=$HAVE_I1, lane=$LANE)."

  {
    printf '%s\n' "$HEADER"
    printf '@RG\tID:%s\tSM:%s\tPU:%s\tPL:ILLUMINA\n' "$RGID" "$RUN" "$RGID"
    printf '@CO\t10x_bam_to_fastq:R1(CR:CY,%s:%s)\n' "$UMI_TAG" "$UMI_QUAL"
    printf '@CO\t10x_bam_to_fastq:R2(SEQ:QUAL)\n'
    if [[ $HAVE_I1 == yes ]]
    then
      printf '@CO\t10x_bam_to_fastq:I1(BC:QT)\n'
    fi
  } > "$RUN.patched.header.sam"

  ## reheader only rewrites the header block and copies the compressed records through;
  ## addreplacerg then stamps RG:Z onto each record. Piping keeps it to one output file.
  $CMD samtools reheader "$RUN.patched.header.sam" "$RUN.bam" \
    | $CMD samtools addreplacerg -R "$RGID" -m overwrite_all --threads "$CPUS" -o "$RUN.patched.bam" -

  rm -f "$RUN.patched.header.sam"
}

function bam2fastq {
  local RUN=$1
  local CPUS=$2
  local CMD=${3:-""}
  local INPUT=$RUN.bam
  local PATCHED=no

  ## this has to be 10x bamtofastq, ideally the latest version
  local HEADER
  HEADER=$($CMD samtools view -H "$RUN.bam")

  if ! header_describes_10x_layout "$HEADER"
  then
    patch_legacy_10x_bam "$RUN" "$CPUS" "$CMD" "$HEADER"
    INPUT=$RUN.patched.bam
    PATCHED=yes
  fi

  $CMD bamtofastq --nthreads $CPUS $INPUT $RUN

  ## The rewritten BAM is as large as the original and nothing downstream reads it.
  if [[ $PATCHED == yes ]]
  then
    rm -f "$INPUT"
  fi

  ## if the chemistry version is v1, the reads would be split into R1 (biological), R2 (barcode, 14 bp) and R3 (UMI, 10 bp).
  ## as we won't be able to process those, we need to fix it

  if [[ `find $RUN/* | grep "_R3_.*fastq.gz"` != "" ]]
  then
    >&2 echo "WARNING: bamtofastq generated read R3! This doesn't usually happen for GEX samples, except for v1 chemistry."
    R1=`find $RUN/* | grep "_R1_.*fastq.gz"`
    R2=`find $RUN/* | grep "_R2_.*fastq.gz"`
    R3=`find $RUN/* | grep "_R3_.*fastq.gz"`
    R1F=`echo $R1 | cut -d' ' -f1`
    R2F=`echo $R2 | cut -d' ' -f1`
    R3F=`echo $R3 | cut -d' ' -f1`
    L1=`zcat $R1F | head -2 | tail -1 | awk '{print length($0)}'`
    L2=`zcat $R2F | head -2 | tail -1 | awk '{print length($0)}'`
    L3=`zcat $R3F | head -2 | tail -1 | awk '{print length($0)}'`
    if (( $L1 > 50 && $L2 == 14 && $L3 == 10))
    then
      >&2 echo "WARNING: v1 chemistry confirmed (read length: R1:$L1, R2:$L2, R3:$L3)! Will concatenate BC+UMI, and move BC+UMI to R1, and biological read to R2.."
      ## inflate the reads - it's more robust than process subs..
      for r2 in $R2
      do
        r3=`echo $r2 | sed "s/_R2_/_R3_/"`
        gzip -d $r2 &
        gzip -d $r3 &
      done
      wait

      ## paste BC+UMI, gzip it
      for r2 in $R2
      do
        r2fq=`echo $r2 | sed "s/\.gz$//"`
        r3fq=`echo $r2 | sed "s/_R2_/_R3_/" | sed "s/\.gz$//"`
        paste $r2fq $r3fq | awk -F '\t' '{if (NR%2==1) {print $1} else {print $1$2}}' | pigz > $r2.fixed &
      done
      wait

      ## move reads around
      for r2 in $R2
      do
        r1=`echo $r2 | sed "s/_R2_/_R1_/"`
        mv $r1 $r2
        mv $r2.fixed $r1
      done

      ## now remove all of the ungzipped fastqs
      rm `find $RUN/* | grep "fastq$"`
    fi
  fi
}

function main () {
  local PARSED=$1
  local CPUS=16
  local RUN=`grep -w "BAM$" $PARSED | cut -f1 | head -$LSB_JOBINDEX | tail -1`
  local SIF="/nfs/cellgeni/singularity/images/reprocess_10x.sif"
  local CMD="singularity run --bind /nfs,/lustre $SIF"

  bam2fastq $RUN $CPUS $CMD
}


if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
