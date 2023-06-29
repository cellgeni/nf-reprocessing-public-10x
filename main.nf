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
    a comma-separated list as a second argument in the sample file  
    (c) Alexander Predeus, Sanger Institute, 2021-2023"
    """.stripIndent()
    exit 1
}

process email_startup {
  
  shell:
  '''
  contents=`cat !{params.samplefile}`
  sendmail "!{params.sangerID}@sanger.ac.uk" <<EOF
  Subject: Launched pipeline
  From: noreply-cellgeni-pipeline@sanger.ac.uk
  Hi there, you've launched Cellular Genetics Informatics' Reprocessing pipeline.
  Your parameters are:
  Samplefile: !{params.samplefile}
  Run STARsolo: !{params.run_starsolo}
  Keep BAMs: !{params.keep_bams}
  Sorting BAM memory (in bytes): !{params.sort_bam_mem}
  
  Your sample file looks like:
  $contents
  Thanks,
  Cellular Genetics Informatics
  EOF
  '''
}

process step1 {
 
  publishDir "/lustre/scratch127/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/reprocessing-results/metadata", mode: "copy"

  input:
  val(sample)
 
  output:
  //Want to copy all metadata to output but need specific files for later processes
  path("*")
  tuple env(SERIES), path("*/*.sample_x_run.tsv"), path("*/*.parsed.tsv"), emit: series_metadata
  path("*/*.parsed.tsv"), emit: series_samples_urls_tsv
  path("*/*.sample.list"), emit: series_sample_list
  path("*/*.run.list"), emit: series_run_list
  path("*/*.sample_x_run.tsv"), emit: series_sample_run_tsv

  shell:
  '''
  SERIES=`echo !{sample} | cut -f 1 -d " "`
  SUBSET=`echo !{sample} | cut -f 2 -d " "`

  ### Step 1: set everything up, collect metadata using ffq/curl from ENA, parse the outputs using custom scripts
  if [[ -d "${SERIES}" ]]; then
    >&2 echo "ERROR: Cannot create directory ${SERIES} because it already exists!"
    exit
  else
    mkdir ${SERIES}
  fi

  cd ${SERIES}
 
  #If there is not sample list path present in samplefile then both variables are set to series ID, this is fixed below
  if [[ "${SERIES}" == "${SUBSET}" ]]; then
    SUBSET=""
  else
    echo $SUBSET | tr "," "\n" > subset.txt
  fi

  if [[ "${SUBSET}" != "" ]]; then
    >&2 echo "WARNING: Using ${SUBSET} to only process select samples!"
    if [[ `grep "^GSM" subset.txt` == "" && `grep "^SRS" subset.txt` == "" && `grep "^ERS" subset.txt` == "" ]]; then
      >&2 echo "ERROR: The subset file ${SUBSET} can only contain GSM, SRS, or ERS IDs!"
      exit 1
    fi
  fi

  cp !{projectDir}/bin/curl_ena_metadata.sh .
  
  if [[ -f subset.txt ]]; then
    !{projectDir}/bin/collect_metadata.sh ${SERIES} subset.txt
  else
    !{projectDir}/bin/collect_metadata.sh ${SERIES} ${SUBSET}
  fi
  ## finally, classify each run into 3 major types:
  ## 1) we have useable 10x paired-end files; 2) we need to get them from 10x BAM; 3) we need to get them from SRA
  ## simultaneously, '$SERIES.urls.list' is generated listing all things that need to be downloaded
  !{projectDir}/bin/parse_ena_metadata.sh $SERIES | sed "s/^/${SERIES}\t/g" > $SERIES.parsed.tsv
  
  ## do not want curl_ena_mtadata.sh script in output so once used, remove it
  rm curl_ena_metadata.sh

  #.command.env is written inside CWD, not work directory set up by nextflow
  #see issue: https://github.com/nextflow-io/nextflow/issues/2812
  #solution is to change back up a directory
  cd ..
  '''
}

process step2 {

  input:
  env(metadata)

  output:
  env(SERIES)
  env(SAMPLE)
  path("done_wget/*")

  shell:
  '''
  SERIES=`echo ${metadata} | cut -f 1 -d " "`
  SAMPLE=`echo ${metadata} | cut -f 2 -d " "`
  URL_INPUT=`echo ${metadata} | cut -f 5 -d " "`
  FILETYPE=`echo ${metadata} | cut -f 6 -d " "`
  # download all the necessary raw files using 'transfer' queue on Farm. Did we tell you this whole thing is Sanger-specific?
  COUNT=1

  # split up fastq files if there are more than one into an array
  if grep -q ";" <<< "${URL_INPUT}"; then
    arrURL=(${URL_INPUT//;/ }) 
  else
    arrURL=(${URL_INPUT})
  fi

  for URL in ${arrURL[@]}; do
    echo "Iteration 1: downloading the files..."
    echo "--------------------------------------------------------" 

    !{projectDir}/bin/wget_restart.sh ${URL}

    ## they did finish! let's run the cleanup  
    echo "Cleanup 1: running the cleanup script..."
    echo "--------------------------------------------------------" 

    !{projectDir}/bin/cleanup_wget_downloads.sh ${URL}

    ## let us repeat until no URLs are left in missing_URLs.list
    while [[ -f missing_URLs.list ]]
    do
      COUNT=$((COUNT+1))
      echo "Iteration $COUNT: downloading the files..." 
      echo "--------------------------------------------------------" 

      missing_URL=`cat missing_URLs.list` #should be a list of 1
      !{projectDir}/bin/wget_restart.sh ${URL}

    ## they did finish! let's run the cleanup again
      echo "Cleanup $COUNT: running the cleanup script..."
      echo "--------------------------------------------------------" 
      !{projectDir}/bin/cleanup_wget_downloads.sh ${URL}
    done

    if [ "${FILETYPE}" == "BAM" ]; then
      mv done_wget/*.bam* "done_wget/${SAMPLE}.bam"
    elif [ "${FILETYPE}" == "GZ1" ] || [ [ "${FILETYPE}" == "GZ2" ]; then
      echo "need to do move fastqs to correct extension too but need to see how fastqs are passed to metadata"
      echo "that could change how this process works for fastqs"
    fi
  done

  echo "FILE DOWNLOAD: ALL DONE!"
  '''
}

process step3 {

  input:
  env(SERIES)
  env(SAMPLE)
  path(infiles)

  output:
  tuple env(SERIES), path("*/*.fastq.gz")

  shell:
  '''
  filetype=`echo !{infiles} | cut -f 1 -d " "`
  filename=`basename ${filetype}`
  extension="${filename##*.}"
  if [[ "${extension}" == "bam" ]]; then
  ## this has to be 10x bamtofastq, ideally the latest version
    bamtofastq --nthreads 16 ${filetype} output
    mv output/*/*.fastq.gz .
    rm -rf output
    for i in *.fastq.gz; do mv $i "${SAMPLE}${i#bamtofastq}"; done
  fi
  ## if extension = basename then no extension 
  if [[ "${extension}" == "${filename}" ]]; then #assume SRA if no extension
    !{projectDir}/bin/sra_to_10x_fastq_gz.sh ${filetype}
  fi
  ## if filename contains "fastq" or "fq" then treat it as a fastq
  if [[ "${filename}" == *"fastq"* ]] || [[ "${filename}" == *"fq"* ]]; then
    for fq in *.f*q*; do
      !{projectDir}/bin/sorting_fastqs.sh "${SAMPLE}" "${fq}"
    done
  fi
  # move to SERIES directory
  mkdir "${SERIES}"
  mv *.gz "${SERIES}"
  # compress any non-compressed fastqs
  find "${SERIES}"/* | while read i; do if [[ $i != *".gz" ]]; then gzip $i; fi; done
  '''
}

process step4 {

  publishDir "/lustre/scratch127/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/reprocessing-results", mode: "copy"

  input:
  path(sample_list)
  path(run_list)
  path(sample_run_list)
  tuple env(SERIES), path(in_fqs)


  output:
  env(SERIES), emit: id //This is only done to ensure email finish happens after workflow completes
  tuple env(SERIES), path("fastqs/*/*"), emit: org_fq
  
  shell:
  '''
  mkdir -p "fastqs/${SERIES}"
  
  SAMPLES=`cat "!{sample_list}"`
  RUNS=`cat "!{run_list}"`
  SAMPLE_RUNS="!{sample_run_list}"

  for i in $RUNS
  do
    if [ ! -z "$(ls "${i}_*.fastq.gz" 2>/dev/null)" ]
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
      if [ ! "$(ls "${j}_*.fastq.gz" 2>/dev/null)" ]
      then
        mv ${j}_*.fastq.gz $i
      fi
    done
    mv $i "fastqs/${SERIES}"
    >&2 echo "Moving directory $i to ./fastqs/${SERIES}"
  done

  echo "REORGANISE FASTQS: ALL DONE!"
  '''
}

process step5 {

  when:
  params.run_starsolo


  input:
  tuple env(SERIES), path(fastq_dir), path(series_runs), path(series_tsv)

  output:
  path("*_starsolo"), emit: ss
  path("*.qc.txt"), emit: qc
  
  shell:
  '''
  SAMPLE=`echo !{fastq_dir}`
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_10x_auto.sh "!{fastq_dir}" "${SERIES}" "${SAMPLE}" "!{series_runs}" "!{series_tsv}" "true" "!{params.sort_bam_mem}"
  else
    !{projectDir}/bin/starsolo_10x_auto.sh "!{fastq_dir}" "${SERIES}" "${SAMPLE}" "!{series_runs}" "!{series_tsv}" "false" "!{params.sort_bam_mem}"
  fi
  !{projectDir}/bin/solo_QC.sh "${SAMPLE}_starsolo" > "${SAMPLE}.qc.txt"
  '''
}

process step6 {

  when:
  params.run_starsolo

  input:
  path(starsolo_dirs)
  path(qc_files)

  output:
  path('starsolo_results/qc.tsv'), emit: qc //This is only done to ensure email finish happens after workflow completes
  path('starsolo_results/*'), emit: results

  shell:
  '''
  #create singular file with qc of every sample
  echo -e "Sample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u" > qc.tsv
  tail -n 1 *qc.txt | grep starsolo >> qc.tsv
  mkdir starsolo_results
  mv *_starsolo starsolo_results
  mv qc.tsv starsolo_results
  '''
}

process email_finish {

  input:
  val(id)

  shell:
  '''
  sendmail "!{params.sangerID}@sanger.ac.uk" <<EOF
  Subject: Finished pipeline
  From: noreply-cellgeni-pipeline@sanger.ac.uk
  Hi there, your run of Cellular Genetics Informatics' reprocessing pipeline is complete.
  Results are available here: "/lustre/scratch127/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/reprocessing-results"
  The results will be deleted in a week so please copy your data to a sensible location, i.e.:
  cp -r "/lustre/scratch127/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/reprocessing-results" /path/to/sensible/location
  Thanks,
  Cellular Genetics Informatics
  EOF
  '''
}

workflow { 
  email_startup()
  ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
  ch_sample_list | flatMap{ it.readLines() } | step1 
  //parallelised downloading urls, ensure all downloads are complete before moving on, ensures all fastqs are generated before running step4 and step5
  //TO DO: change pipeline so STARsolo runs when each run's fastqs are available, not all the fastqs
  //below step groups fastqs series id, then creates tuple of series id and all fastqs which are passed to step4
  step1.out.series_samples_urls_tsv | splitText | step2 | step3 | groupTuple | map ({ it -> [it.first(), it.tail().flatten()] }) | set { ch4 }
  step4(step1.out.series_sample_list, step1.out.series_run_list, step1.out.series_sample_run_tsv, ch4)
  //below adds the series_samples_urls_tsv and series_sample_run_tsv to each set of fastqs before transposing so each tuple is its own channel
  if (params.run_starsolo == false) {
    email_finish(step4.out.id)
  }
  step4.out.org_fq | join( step1.out.series_metadata ) | transpose | step5
  step5.out.ss | collect | set { step5_ss }
  step5.out.qc | collect | set { step5_qc }
  step6(step5_ss, step5_qc) 
  step6.out.results | flatten |  subscribe { it -> itname = it.getName(); it.copyTo("/lustre/scratch127/cellgen/cellgeni/tickets/nextflow-tower-results/${params.sangerID}/${params.timestamp}/reprocessing-results/starsolo/${itname}") }
  //Ternary operator to ensure finishing email is sent after final process depending on whether starsolo is ran or not
  params.run_starsolo == true ? email_finish(step6.out.qc) : email_finish(step4.out.id 
}
