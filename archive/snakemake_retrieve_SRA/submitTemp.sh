#!/bin/bash
#$ -N sra
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra
#$ -pe smp 2
#$ -l mem=5G,time=12::
#
#######################################################################
# Load the conda environment for running snakemake on the master job
echo "Loading conda"
module load conda
echo "Activating snakemake env"
source activate snakemake
#######################################################################
# The target rule for the pipeline
TARGETRULE="cosma_IDs"
#######################################################################
# The logfile for jobs
LOGFILE="/ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra/status.txt"
#######################################################################
# Set up the desired variables for running the job
echo "Setting up variables"
SNAKEFILE="/ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra/retrieve_SRA.snakefile"
WORKDIR="/ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra"
CONDADIR="conda"
CLUSTCONFIG="/ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra/clustConfig.json"
CORES="999"
LOCALCORES="1"
JOBRESTARTS="5"
LATENCYWAIT="120"
#######################################################################
# Initiate the pipeline
echo "Spooling the snakemake pipeline"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
--local-cores $LOCALCORES \
$TARGETRULE \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -terse -j y -o /ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra/temp.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/cosma_sra -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--max-jobs-per-second 1000
echo "Finished"
