#!/bin/bash
#$ -N wgs_lncap_fork1
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
#$ -pe smp 4
#$ -l mem=5G,time=24::
#
# Load the conda environment for running snakemake on the master job
echo "Loading conda"
module load conda
echo "Activating snakemake env"
source activate snakemake
# Set up the desired variables for running the job
echo "Setting up variables"
#TARGETRULE="footprinting_wt01_pwm95"
SNAKEFILE="/ifs/scratch/c2b2/ac_lab/jk3755/atac/ATAC.snakefile"
WORKDIR="/ifs/scratch/c2b2/ac_lab/jk3755/atac"
CONDADIR="conda"
CLUSTCONFIG="/ifs/scratch/c2b2/ac_lab/jk3755/atac/resources/config/clustConfig.json"
CORES="100"
LOCALCORES="4"
JOBRESTARTS="5"
LATENCYWAIT="120"
# Echo variable settings
echo "Current cluster settings for snakemake run:"
echo "Target rule: $TARGETRULE"
echo "Snakefile: $SNAKEFILE"
echo "Working directory: $WORKDIR"
echo "Conda environment directory: $CONDADIR"
echo "Cluster configuration file: $CLUSTCONFIG"
echo "Provided cluster cores: $CORES"
echo "Provided local cores: $LOCALCORES"
echo "Job restart attempts: $JOBRESTARTS"
echo "Filesystem latency wait: $LATENCYWAIT"
#### Run the pipeline, single logfile
echo "Spooling the snakemake pipeline"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
--local-cores $LOCALCORES \
wgs_lncap_fork1 \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -terse -j y -o /ifs/scratch/c2b2/ac_lab/jk3755/atac/wgs_lncap_fork1.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--max-jobs-per-second 1000
