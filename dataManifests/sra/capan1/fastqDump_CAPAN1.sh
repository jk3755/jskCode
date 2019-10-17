#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## CAPAN1 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5437544[accn]
bin/fastq-dump --outdir capan1/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8639189

## CAPAN1 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432692[accn]
bin/fastq-dump --outdir capan1/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633979	

## CAPAN1 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414506[accn]
bin/fastq-dump --outdir capan1/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615247