#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## KP4 WGS: https://www.ncbi.nlm.nih.gov/sra/ERX2438302[accn]
bin/fastq-dump --outdir kp4/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip ERR2397499

## KP4 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432751[accn]
bin/fastq-dump --outdir kp4/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633920	

## KP4 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414975[accn]
bin/fastq-dump --outdir kp4/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8616075