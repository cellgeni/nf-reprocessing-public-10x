#!/bin/bash -e

SAMPLE=$1
fq=$2

regex="^(.+)_(S[0-9]+)_([L,R][0-9]+)_([R,I][0-9]+)_([0-9]+)\.fastq(\.gz)?$"
prefix="${fq%%.*}"
# check if fastq naming convention matches cellranger required naming
if [[ "${fq}" =~ $regex ]]; then
  # if so then remove the prefix before and replace it with the sample ID
  mv $fq "${SAMPLE}${fq#${BASH_REMATCH[1]}}"
else
  # if not then check if the fastq naming convention is "_1", "_2" or "_R1", "_R2", "_I1", "_I2"
  regex2="_([R,I]?[1-2])"
  if [[ "${fq}" =~ $regex2 ]]; then
    if [[ $(echo ${BASH_REMATCH[1]} | grep "[R,I][1,2]") ]]; then
      mv $fq "${SAMPLE}_S1_L001_${BASH_REMATCH[1]}_001${fq#${prefix}}"
    elif [[ $(echo ${BASH_REMATCH[1]} | grep "[1,2]") ]]; then
      mv $fq "${SAMPLE}_S1_L001_R${BASH_REMATCH[1]}_001${fq#${prefix}}"
    fi
  fi
fi
