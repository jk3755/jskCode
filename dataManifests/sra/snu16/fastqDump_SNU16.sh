#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## SNU16 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5466612[accn]
bin/fastq-dump --outdir snu16/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8670768

## SNU16 HC: https://www.ncbi.nlm.nih.gov/sra/SRX5454943[accn]
bin/fastq-dump --outdir snu16/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8657217

## SNU16 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432497[accn]
bin/fastq-dump --outdir snu16/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633575

## SNU16 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5415116[accn]
bin/fastq-dump --outdir snu16/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615934