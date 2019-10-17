#!/bin/bash
##
## Use parallel fastq-dump: https://github.com/rvalieris/parallel-fastq-dump
##
## H508 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5466690[accn]
parallel-fastq-dump --sra-id SRR8670690 --threads 20 --outdir h508 --split-files --gzip

## SNU16 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5466612[accn]
parallel-fastq-dump --sra-id SRR8670768 --threads 20 --outdir snu16/

## CAPAN1 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5437544[accn]
parallel-fastq-dump --sra-id SRR8639189 --threads 20 --outdir capan1/

## LNCAP WGS: https://www.ncbi.nlm.nih.gov/sra/SRX2541290[accn]
parallel-fastq-dump --sra-id SRR5233717 --threads 20 --outdir lncap/
parallel-fastq-dump --sra-id SRR5259501 --threads 20 --outdir lncap/
parallel-fastq-dump --sra-id SRR5259502 --threads 20 --outdir lncap/
parallel-fastq-dump --sra-id SRR5259503 --threads 20 --outdir lncap/

## PANC1 WGS: https://www.ncbi.nlm.nih.gov/sra/?term=panc1
parallel-fastq-dump --sra-id SRR8670733 --threads 20 --outdir panc1/