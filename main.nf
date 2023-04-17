#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def errorMessage() {
    log.info"""
    ==================
    reprocessing error
    ==================
    USAGE: nextflow run main.nf --samplefile /path/to/samplefile
    The samplefile should be a tab separated file containing a series ID followed by an optional path to a list of samples within that series.
    This script should take a series ID (GSE from GEO, E-MTAB from ArrayExpress, or PRJ* from bioProjects/SRA),
    pull/parse all relevant metadata, download/reformat all raw read files, and return fully reprocessed count 
    (I won't promise you it will work every time, but I promise you it will try!)
    In case you don't want all the samples from the series processed, provide a list of samples you do want as 
    (c) Alexander Predeus, Sanger Institute, 2021-2023"
    """.stripIndent()
    exit 1
}

process step1 {
 
  publishDir "${params.outdir}/metadata", mode: "copy"

  input:
  val(sample)
  
  output:
  env(SERIES), emit: series_id
  path("*/*.parsed.tsv"), emit: series_samples_urls_tsv
  path("*/*.sample.list"), emit: series_sample_list
  path("*/*.run.list"), emit: series_run_list
  path("*/*.sample_x_run.tsv"), emit: series_sample_run_tsv

  shell:
  '''
  SERIES=`echo !{sample} | cut -f 1 -d " "`
  SUBSET=`echo !{sample} | cut -f 2 -d " "`
 
  #If there is not sample list path present in samplefile then both variables are set to series ID, this is fixed below
  if [[ "${SERIES}" == "${SUBSET}" ]]; then
    SUBSET=""
  fi

  ### Step 1: set everything up, collect metadata using ffq/curl from ENA, parse the outputs using custom scripts

  if [[ -d "${SERIES}" ]]; then
    >&2 echo "ERROR: Cannot create directory ${SERIES} because it already exists!"
    exit
  else
    mkdir ${SERIES}
  fi

  if [[ "${SUBSET}" != "" ]]; then
    >&2 echo "WARNING: Using file ${SUBSET} to only process select samples!"
    if [[ `grep "^GSM" ${SUBSET}` == "" && `grep "^SRS" ${SUBSET}` == "" && `grep "^ERS" ${SUBSET}` == "" ]]; then
      >&2 echo "ERROR: The subset file ${SUBSET} can only contain GSM, SRS, or ERS IDs!"
      exit 1
    fi
  fi

  cd ${SERIES}
  cp !{projectDir}/bin/curl_ena_metadata.sh .
  
  !{projectDir}/bin/collect_metadata.sh ${SERIES} ${SUBSET}
  ## finally, classify each run into 3 major types:
  ## 1) we have useable 10x paired-end files; 2) we need to get them from 10x BAM; 3) we need to get them from SRA
  ## simultaneously, '$SERIES.urls.list' is generated listing all things that need to be downloaded
  !{projectDir}/bin/parse_ena_metadata.sh $SERIES > $SERIES.parsed.tsv
  
  #.command.env is written inside CWD, not work directory set up by nextflow
  #see issue: https://github.com/nextflow-io/nextflow/issues/2812
  #solution is to change back up a directory
  cd ..
  '''
}

process step2 {

  input:
  tuple env(metadata), env(SERIES)

  output:
  env(SAMPLE)
  path("done_wget/*")

  shell:
  '''
  SAMPLE=`echo ${metadata} | cut -f 1 -d " "`
  url_path=`echo ${metadata} | cut -f 4 -d " "`
  ## download all the necessary raw files using 'transfer' queue on Farm. Did we tell you this whole thing is Sanger-specific?
  COUNT=1
  ## since this is a part of reprocessing script, we assume all the necessary sub-scripts have been copied already

  echo "Iteration 1: downloading the files..."
  echo "--------------------------------------------------------" 

  !{projectDir}/bin/wget_restart.sh ${url_path}

  ## they did finish! let's run the cleanup  
  echo "Cleanup 1: running the cleanup script..."
  echo "--------------------------------------------------------" 
  !{projectDir}/bin/cleanup_wget_downloads.sh ${url_path}

  ## let us repeat until no URLs are left in missing_URLs.list
  while [[ -f missing_URL.list ]]
  do
    COUNT=$((COUNT+1))
    echo "Iteration $COUNT: downloading the files..." 
    echo "--------------------------------------------------------" 
  
    missing_URL=`cat missing_URLs.list` #should be a list of 1
    !{projectDir}/bin/wget_restart.sh ${missing_URL}
   
  ## they did finish! let's run the cleanup again
    echo "Cleanup $COUNT: running the cleanup script..."
    echo "--------------------------------------------------------" 
    !{projectDir}/bin/cleanup_wget_downloads.sh ${url_path}
  done

  echo "FILE DOWNLOAD: ALL DONE!" 
  '''
}

process step3 {

  input:
  env(SAMPLE)
  path(BAM_or_SRA)

  output:
  path("*.fastq.gz")

  shell:
  '''
  filename=`basename !{BAM_or_SRA}`
  extension="${filename##*.}"
  if [[ “${extension}” == “bam” ]]; then
  ## this has to be 10x bamtofastq, ideally the latest version
    bamtofastq --nthreads 16 !{BAM_or_SRA} output
    mv output/*/*.fastq.gz .
    rm -rf output
    for i in *.fastq.gz; do mv $i “${SAMPLE}${i#bamtofastq}“; done
  fi
  ## if extension = basename then no extension 
  if [[ “${extension}” == “${filename}” ]]; then #assume SRA if no extension
    !{projectDir}/bin/sra_to_10x_fastq_gz.sh !{BAM_or_SRA}
  fi
  '''
}

process step4 {

  publishDir "${params.outdir}", mode: "copy"


  input:
  env(SERIES)
  path(fastqs)

  output:
  path("organised_fastqs"), emit: fastqs

  shell:
  '''
  mkdir organised_fastqs

  SAMPLES=`cat "!{params.outdir}/metadata/${SERIES}/${SERIES}.sample.list"`
  RUNS=`cat "!{params.outdir}/metadata/${SERIES}/${SERIES}.run.list"`
  SAMPLE_RUNS="!{params.outdir}/metadata/${SERIES}/${SERIES}.sample_x_run.tsv"
  for i in $RUNS
  do
    if [[ ! -s "${i}_1.fastq.gz" || ! -s "${i}_2.fastq.gz" ]] 
    then
      >&2 echo "WARNING: Run $i does not seem to have two fastq files (or a bamtofastq output directory) associated with it! Please investigate."
    fi
  done

  for i in $SAMPLES
  do
    >&2 echo "Moving the files for sample $i:"
    mkdir $i
    SMPRUNS=`grep $i ${SAMPLE_RUNS} | awk '{print $2}' | tr ',' '\n'`
    for j in $SMPRUNS
    do
      >&2 echo "==> Run $j belongs to sample $i, moving to directory $i.."
      if [[ -s "${j}_1.fastq.gz" ]]
      then
        mv ${j}_?.fastq.gz $i
      fi
    done
    mv $i organised_fastqs
    >&2 echo "Moving directory $i to /fastqs.."
  done

  echo "REORGANISE FASTQS: ALL DONE!"
  '''
}

process step5 {

  publishDir "${params.outdir}/starsolo_results", mode: "copy", pattern: "*_starsolo"

  input:
  tuple env(SAMPLE), path(organised_fastqs), env(SERIES)

  output:
  path("*_starsolo")
  path("*.qc.txt"), emit: qc
  
  
  shell:
  '''
  trimmed_SAMPLE=`echo $SAMPLE | tr -d '\n'`
  !{projectDir}/bin/starsolo_10x_auto.sh "!{organised_fastqs}" $SERIES ${trimmed_SAMPLE} "!{params.outdir}/metadata/${SERIES}" 
  !{projectDir}/bin/solo_QC.sh "${trimmed_SAMPLE}_starsolo" > "${trimmed_SAMPLE}.qc.txt"
  '''
}

process step6 {

  publishDir "${params.outdir}/starsolo_results", mode: "copy"

  input:
  path(qc_files)

  output:
  path('qc.tsv')

  shell:
  '''
  #create singular file with qc of every sample
  echo -e "Sample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u" > qc.tsv
  tail -n 1 *qc.txt | grep starsolo >> qc.tsv
  '''
}

workflow {
  ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
  ch_sample_list | flatMap{ it.readLines() } | step1 
  //parallelised downloading urls, ensure all downloads are complete before moving on, ensures all fastqs are generated before running step4 and step5
  //TO DO: change pipeline so STARsolo runs when each run's fastqs are available, not all the fastqs
  step1.out.series_samples_urls_tsv | splitText | combine(step1.out.series_id) | step2 | step3 | collect | set {ch_fastqs}
  step4(step1.out.series_id, ch_fastqs)
  step1.out.series_sample_list | splitText | combine(step4.out.fastqs) | combine(step1.out.series_id) | step5
  step5.out.qc | collect | step6
}
