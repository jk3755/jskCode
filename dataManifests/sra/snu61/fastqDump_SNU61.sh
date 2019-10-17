#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## SNU-61: WXS https://www.ncbi.nlm.nih.gov/sra/SRX5418100[accn] Accession number SRR8618987
bin/fastq-dump --outdir snu61/wxs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8618987

## SNU-61: RNAseq https://www.ncbi.nlm.nih.gov/sra/SRX5415026[accn] Accession number SRR8616024
bin/fastq-dump --outdir snu61/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8616024

## SNU-61: RRBS https://www.ncbi.nlm.nih.gov/sra/SRX5432783[accn]
bin/fastq-dump --outdir snu61/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633888
