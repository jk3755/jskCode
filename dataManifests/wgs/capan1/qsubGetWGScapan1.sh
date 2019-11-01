#!/bin/bash
#$ -N get.CAPAN1.WGS
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l mem=25G,time=48:0:0
#
# https://www.ncbi.nlm.nih.gov/sra/SRX5437544[accn]
#
echo "Loading conda"
module load conda
echo "Activating parallel-fastq-dump env"
source activate parallel-fastq-dump
#
echo "SRR8639189 prefetch"
prefetch -v SRR8639189 -X 100G
echo "SRR8639189 fastq dump"
fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/capan1/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR8639189/SRR8639189.sra
