#!/bin/bash
#$ -N get.SNU16.WGS
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l mem=25G,time=48:0:0
#
# https://www.ncbi.nlm.nih.gov/sra/SRX5466612[accn]
#
echo "Loading conda"
module load conda
echo "Activating parallel-fastq-dump env"
source activate parallel-fastq-dump
#
echo "SRR8670768 prefetch"
prefetch -v SRR8670768 -X 100G
echo "SRR8670768 fastq dump"
fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/snu16/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR8670768/SRR8670768.sra
