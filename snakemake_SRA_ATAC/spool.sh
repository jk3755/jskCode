#!/bin/bash
#
#
snakemake --snakefile SRA_ATAC.snakefile --cores 20 processID --use-conda --conda-prefix conda --restart-times 10 --latency-wait 120 --rerun-incomplete