#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

## LNCaP WGS: https://www.ncbi.nlm.nih.gov/sra/SRX2541290[accn]
bin/fastq-dump --outdir lncap/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR5259503
bin/fastq-dump --outdir lncap/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR5259502
bin/fastq-dump --outdir lncap/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR5259501
bin/fastq-dump --outdir lncap/wgs --gzip --skip-technical --readids --dumpbase --split-files --clip SRR5233717
