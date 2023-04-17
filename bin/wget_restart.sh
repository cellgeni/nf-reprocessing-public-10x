#!/bin/bash 

URL=$1

TAG=`basename ${URL%%.fastq.gz}` 
wget --retry-connrefused --read-timeout=20 --timeout=15 --tries=0 --continue $URL &> "wget.log"
