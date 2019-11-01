#!/bin/bash
#$ -N seqkitSplitLNCAP
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l mem=20G,time=6:0:0
#
echo "Loading conda"
module load conda
echo "Activating seqkit env"
source activate seqkit
#
echo "Split SRR5233717 in 20 parts"
seqkit split2 -p 20 -f -O lncap/split1 -1 lncap/SRR5233717_1.fastq -2 lncap/SRR5233717_2.fastq
#
echo "Split SRR5259501 in 20 parts"
seqkit split2 -p 20 -f -O lncap/split2 -1 lncap/SRR5259501_1.fastq -2 lncap/SRR5259501_2.fastq
#
echo "Split SRR5259502 in 20 parts"
seqkit split2 -p 20 -f -O lncap/split3 -1 lncap/SRR5259502_1.fastq -2 lncap/SRR5259502_2.fastq
#
# This one doesnt work
#echo "Split SRR5259503 in 20 parts"
#seqkit split2 -p 20 -f -O lncap/split -1 lncap/SRR5259503_1.fastq -2 lncap/SRR5259503_2.fastq