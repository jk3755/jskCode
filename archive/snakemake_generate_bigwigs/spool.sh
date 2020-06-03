#!/bin/bash
#
#
snakemake --snakefile generate_bigwigs.snakefile --cores 20 sample_expander --use-conda --conda-prefix conda --restart-times 5 --latency-wait 120 --rerun-incomplete