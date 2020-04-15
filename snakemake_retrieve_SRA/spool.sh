#!/bin/bash
#
#
snakemake --snakefile retrieve_SRA.snakefile --cores 20 cosma_IDs --restart-times 5 --latency-wait 60 --keep-going