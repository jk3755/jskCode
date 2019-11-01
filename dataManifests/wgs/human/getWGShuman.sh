#!/bin/bash
######################################################
## https://www.ncbi.nlm.nih.gov/sra/?term=SRR1554094
## sra-tools: https://github.com/ncbi/sra-tools
######################################################
## Human WGS using CPT-seq
## Acc. SRR1554094
fastq-dump --split-files -O data SRR1554094