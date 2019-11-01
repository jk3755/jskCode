#!/bin/bash
################################################################################
## https://www.ncbi.nlm.nih.gov/sra/SRX2541290[accn]
## sra-tools: https://github.com/ncbi/sra-tools
################################################################################
## Alternatively, you can use parallel-fastq-dump
## This can significantly reduce download times
## parallel-fastq-dump: https://github.com/rvalieris/parallel-fastq-dump
################################################################################
## LNCaP WGS 
## Four HiSeq runs
## Acc. SRR5233717 SRR5259501 SRR5259502 SRR5259503
fastq-dump --split-files -O data SRR5233717
fastq-dump --split-files -O data SRR5259501
fastq-dump --split-files -O data SRR5259502
fastq-dump --split-files -O data SRR5259503