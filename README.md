# nf-reprocessing-public-10x
Nextflow pipeline of our reprocessing [pipeline](https://github.com/cellgeni/reprocess_public_10x).

## Contents of Repo:
* `main.nf` - the Nextflow pipeline that executes reprocessing.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `examples/samples.list` - example of the expected samplefile containing series IDs and optionally sample IDs.
* `examples/RESUME` - an example run script that executes the pipeline it has 1 hardcoded argument: `/path/to/sample.list` which you will need to change on your installation. 
* `bin` - a directory containing various scripts the pipeline uses to download project data and realign the data with STARsolo.
* `Dockerfile` - a dockerfile to reproduce the environment for step3, step5 has its own container that has its Dockerfile located [here](https://github.com/cellgeni/STARsolo/blob/main/Dockerfile)

## Pipeline Arguments:
* `--samplefile` - The path to the sample file provided to the pipeline which contains one sample ID per line. This sample is assumed to have CRAM files stored on IRODS.
* `--outdir` - The path to where the results will be saved.
* `--run_starsolo` - Tells pipeline whether to realign data with STARsolo or not
* `--keep_bams` - Tells the pipeline whether to generate BAM files (default false means do not generate).
* `--sort_bam_mem` - Input memory (IN BYTES) for starsolo to use for sorting BAM files if BAM files are kept (default 60GB = 60000000000B).
