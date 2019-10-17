#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## PANC0403 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432608[accn]
bin/fastq-dump --outdir panc0403/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8634063	

## PANC0403 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414635[accn]
bin/fastq-dump --outdir panc0403/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615765