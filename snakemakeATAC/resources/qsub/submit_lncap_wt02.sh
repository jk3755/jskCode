#!/bin/bash
#$ -N lncap_wt02
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
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
TARGETRULE="lncap_wt02"
#######################################################################
# The logfile for jobs
LOGFILE="/ifs/scratch/c2b2/ac_lab/jk3755/atac/lncap_wt02.txt"
#######################################################################
# Set up the qsub command string
echo "Setting up qsub command"
Q1="'qsub -terse -j y -o "
Q2="-pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V'"
QSUBCOMMAND=${Q1}${LOGFILE}${Q2}
echo "Qsub command:"
echo $QSUBCOMMAND
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
$TARGETRULE \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -terse -j y -o /ifs/scratch/c2b2/ac_lab/jk3755/atac/lncap_wt02.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--max-jobs-per-second 1000
echo "Finished"
