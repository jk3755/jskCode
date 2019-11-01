#!/bin/bash
#$ -N getLNCaPWGS
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l mem=25G,time=48:0:0
#
echo "Loading conda"
module load conda
echo "Activating sratools env"
source activate sratools
#
# echo "SRR5233717 prefetch"
# prefetch -v SRR5233717 -X 100G
# echo "SRR5233717 fastq dump"
# fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/lncap/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR5233717/SRR5233717.sra
#
# echo "SRR5259501 prefetch"
# prefetch -v SRR5259501 -X 100G
# echo "SRR5259501 fastq dump"
# fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/lncap/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR5259501/SRR5259501.sra
#
# echo "SRR5259502 prefetch"
# prefetch -v SRR5259502 -X 100G
# echo "SRR5259502 fastq dump"
# fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/lncap/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR5259502/SRR5259502.sra
#
# This one doesn't work
# echo "SRR5259503 prefetch"
# prefetch -v SRR5259503 -X 100G
# echo "SRR5259503 fastq dump"
# fastq-dump --outdir /ifs/scratch/c2b2/ac_lab/jk3755/sra/lncap/ --split-files /ifs/scratch/c2b2/ac_lab/jk3755/sra/SRR5259503/SRR5259503.sra
