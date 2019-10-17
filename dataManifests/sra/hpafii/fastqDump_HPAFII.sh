#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## HPAFII RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432464[accn]
bin/fastq-dump --outdir hpafii/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633608	

## HPAFII RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5415193[accn]
bin/fastq-dump --outdir hpafii/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8616209