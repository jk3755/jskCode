#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## MDST8 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414820[accn]
bin/fastq-dump --outdir mdst8/rna --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615580

## MDST8 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5432289[accn]
bin/fastq-dump --outdir mdst8/rrbs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633783