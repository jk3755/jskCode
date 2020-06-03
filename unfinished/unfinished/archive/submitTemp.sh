#!/bin/bash
#$ -N temp
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
#$ -pe smp 4
#$ -l mem=4G,time=12::
#
#######################################################################
# Load the conda environment for running snakemake on the master job
module load conda
source activate snakemake
#######################################################################
# Set up the desired variables for running the job
echo "Setting up variables"
SNAKEFILE="/ifs/scratch/c2b2/ac_lab/jk3755/atac/ATAC.snakefile"
WORKDIR="/ifs/scratch/c2b2/ac_lab/jk3755/atac"
CONDADIR="conda"
CLUSTCONFIG="/ifs/scratch/c2b2/ac_lab/jk3755/atac/resources/config/clustConfig.json"
CORES="999"
LOCALCORES="4"
JOBRESTARTS="5"
LATENCYWAIT="120"
#######################################################################
# Initiate the pipeline
echo "Spooling the snakemake pipeline"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
--local-cores $LOCALCORES \
ciccia_ap1 \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -terse -j y -o /ifs/scratch/c2b2/ac_lab/jk3755/atac/temp.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--max-jobs-per-second 1000
echo "Finished"