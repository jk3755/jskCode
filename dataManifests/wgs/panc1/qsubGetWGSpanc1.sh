#!/bin/bash
#$ -N get.PANC1.WGS
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l mem=25G,time=48:0:0
#
# https://www.ncbi.nlm.nih.gov/sra/SRX5466647[accn]
#
echo "Loading conda"
module load conda
echo "Activating parallel-fastq-dump env"
source activate parallel-fastq-dump
#
echo "SRR8670733 prefetch"
prefetch -v SRR8670733 -X 100G
echo "SRR8670733 fastq dump"
fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/panc1/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR8670733/SRR8670733.sra
