#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## PANC1 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5466647[accn]
bin/fastq-dump --outdir panc1/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8670733

## PANC1 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432767[accn]
bin/fastq-dump --outdir panc1/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633904	

## PANC1 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414715[accn]
bin/fastq-dump --outdir panc1/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615685