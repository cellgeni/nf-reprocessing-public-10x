#!/bin/bash 

if (( $# == 0 ))
then
  echo "USAGE: ./cleanup_wget_downloads.sh <URL>"
  echo "Please provide the URL used for original download!"
  exit 1
fi

URL=$1
TAG=`basename ${URL%%.fastq.gz}` ## this should work for both SRA (SRR1245 or SRR1245.1 formats), and *_1.fastq.gz/*_2.fastq.gz downloads 

if [[ -d done_wget ]]
then
  echo "Found directory named 'done_wget', will use it to move completed downloads!"
else 
  echo "No directory named 'done_wget' found! will create it and use to move completed downloads!"
  mkdir done_wget
fi

if [[ -d logs ]]
then
  echo "Found directory named 'logs', will use it to move completed downloads!"
else
  echo "No directory named 'logs' found! will create it and use to move completed downloads!"
  mkdir logs
fi

if [[ -f missing_URLs.list ]] 
then
  echo "Found file named 'missing_URLS.list'; deleting it.." 
  rm missing_URLs.list
fi

if [[ -f "${TAG}.wget.log" ]]
then 
  SV=`tail "${TAG}.wget.log" | grep -wF saved | wc -l`
  if (( $SV == 1 ))
  then
    echo "Sample $TAG was downloaded successfully, moving files to /done_wget!"
    mv ${TAG}*.wget.log logs
    mv ${TAG}* done_wget
  else 
    AC=$((AC+1)) 
    echo "Sample $TAG was NOT downloaded successufully; removing remaining files and adding URL to 'missing_URLs.list'!" 
    rm ${TAG}*
    echo $TAG >> missing_URLs.list
  fi
else 
  echo "Sample $TAG   --  download was not attempted, probably completed earlier..."
fi
