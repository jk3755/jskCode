#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## PK45H RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432496[accn]
bin/fastq-dump --outdir pk45h/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633576	

## PK45H RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414458[accn]
bin/fastq-dump --outdir pk45h/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615295