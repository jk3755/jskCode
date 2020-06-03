#!/bin/bash
#$ -N cosma_donorA_luf
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
#$ -pe smp 4
#$ -l mem=10G,time=24::
#
#######################################################################
# Load the conda environment for running snakemake on the master job
echo "Loading conda"
module load conda
echo "Activating snakemake env"
source activate snakemake
#######################################################################
# The target rule for the pipeline
TARGETRULE="cosma_donorA_luf"
#######################################################################
# Set up the desired variables for running the job
echo "Setting up variables"
SNAKEFILE="/ifs/scratch/c2b2/ac_lab/jk3755/atac/ATAC.snakefile"
WORKDIR="/ifs/scratch/c2b2/ac_lab/jk3755/atac"
CONDADIR="conda"
CLUSTCONFIG="/ifs/scratch/c2b2/ac_lab/jk3755/atac/resources/config/clustConfig.json"
CORES="150"
LOCALCORES="4"
JOBRESTARTS="3"
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
--cluster "qsub -terse -j y -o /ifs/scratch/c2b2/ac_lab/jk3755/atac/cosma_donorA_luf.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--max-jobs-per-second 1000
echo "Finished"
