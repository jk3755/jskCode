#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

#### H508 #######################################################################################################

## H508 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5466690[accn]
bin/fastq-dump --outdir h508/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8670690

## H508 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5431921[accn]
bin/fastq-dump --outdir h508/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633497	

## H508 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414591[accn]
bin/fastq-dump --outdir h508/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615809
