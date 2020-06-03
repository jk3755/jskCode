#!/bin/bash
#$ -N WGS
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/wgs
#$ -pe smp 1
#$ -l mem=10G,time=48::
#
echo "Loading conda"
module load conda
echo "Activating snakemake env"
source activate snakemake
#
echo "Unlocking working directory"
snakemake \
--snakefile /ifs/scratch/c2b2/ac_lab/jk3755/wgs/WGS.snakefile \
--cores 999 \
UNLOCK \
--cluster-config /ifs/scratch/c2b2/ac_lab/jk3755/wgs/resources/config/clusterConfig.json \
--cluster "qsub -terse -j y -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/wgs -V" \
--use-conda \
--conda-prefix conda \
--restart-times 2 \
--latency-wait 60 \
--rerun-incomplete \
--unlock
#
echo "Spooling the pipeline to cluster"
snakemake \
--snakefile /ifs/scratch/c2b2/ac_lab/jk3755/wgs/WGS.snakefile \
--cores 999 \
human_SRR1554094 \
--cluster-config /ifs/scratch/c2b2/ac_lab/jk3755/wgs/resources/config/clusterConfig.json \
--cluster "qsub -terse -j y -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/wgs -V" \
--use-conda \
--conda-prefix conda \
--restart-times 2 \
--latency-wait 60 \
--rerun-incomplete
