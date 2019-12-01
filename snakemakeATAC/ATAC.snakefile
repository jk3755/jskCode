########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################
configfile: "resources/config/config.yaml"
include: "resources/modules/spoolPreprocessing.snakefile"
include: "resources/modules/spoolFootprinting.snakefile"
include: "resources/modules/spoolFullAnalysis.snakefile"

########################################################################################################################################
#### UNLOCK TARGET #####################################################################################################################
########################################################################################################################################

## Just a dummy rule for quickly unlocking the working directory
## Is this really the only way to do this?
rule UNLOCK:
    output:
        "unlock.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### FULL ANALYSIS #####################################################################################################################
########################################################################################################################################

## This rule determines what is run in the full analysis spooling option
## Also cleans up the intermediate preprocessing data
rule AGGREGATOR_full_analysis:
    input:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.preprocessing",
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.uncorrected"
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.full"
    shell:
        """
        rm -rf {wildcards.path}preprocessing
        touch {output}
        """
        
########################################################################################################################################
#### DIRECTORY STRUCTURE ###############################################################################################################
########################################################################################################################################

rule AGGREGATOR_build_directory_structure:
    input:
        "{path}operations/directories/main_dir.built",
        "{path}operations/directories/operations_dir.built",
        "{path}operations/directories/preprocessing_dir.built",
        "{path}operations/directories/peaks_dir.built",
        "{path}operations/directories/metrics_dir.built",
        "{path}operations/directories/footprints_dir.built" 
    output:
        "{path}operations/directories/all_dirs.built"
    shell:
        "touch {output}"

rule DIR_main:
    output:
        "{path}operations/directories/main_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}correlation
        mkdir -p -v {wildcards.path}figures
        mkdir -p -v {wildcards.path}aligned
        mkdir -p -v {wildcards.path}bigwig
        mkdir -p -v {wildcards.path}footprints
        touch {output}
        """

## Important to have separated dir structure here
rule DIR_operations:
    output:
        "{path}operations/directories/operations_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations/directories
        mkdir -p -v {wildcards.path}operations/aggregators
        touch {output}
        """

## Avoid repeating words in the full paths to intermediate files
rule DIR_preprocessing:
    output:
        "{path}operations/directories/preprocessing_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing/1
        mkdir -p -v {wildcards.path}preprocessing/2
        mkdir -p -v {wildcards.path}preprocessing/3
        mkdir -p -v {wildcards.path}preprocessing/4
        mkdir -p -v {wildcards.path}preprocessing/5
        mkdir -p -v {wildcards.path}preprocessing/6
        mkdir -p -v {wildcards.path}preprocessing/7
        mkdir -p -v {wildcards.path}preprocessing/8
        mkdir -p -v {wildcards.path}preprocessing/9
        mkdir -p -v {wildcards.path}preprocessing/a
        mkdir -p -v {wildcards.path}preprocessing/b
        mkdir -p -v {wildcards.path}preprocessing/c
        mkdir -p -v {wildcards.path}preprocessing/d
        mkdir -p -v {wildcards.path}preprocessing/e
        mkdir -p -v {wildcards.path}preprocessing/f
        mkdir -p -v {wildcards.path}preprocessing/g
        mkdir -p -v {wildcards.path}preprocessing/h
        touch {output}
        """

rule DIR_metrics:
    output:
        "{path}operations/directories/metrics_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}metrics/fastq
        mkdir -p -v {wildcards.path}metrics/myco
        mkdir -p -v {wildcards.path}metrics/align
        mkdir -p -v {wildcards.path}metrics/genomecov
        mkdir -p -v {wildcards.path}metrics/size
        mkdir -p -v {wildcards.path}metrics/fragsizes
        mkdir -p -v {wildcards.path}metrics/peakideogram
        mkdir -p -v {wildcards.path}metrics/duplication
        mkdir -p -v {wildcards.path}metrics/peakanno
        touch {output}
        """

rule DIR_peaks:
    output:
        "{path}operations/directories/peaks_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}peaks/ln {wildcards.path}peaks/gn {wildcards.path}peaks/sm
        touch {output}
        """

rule DIR_footprints:
    output:
        "{path}operations/directories/footprints_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}footprints/ss
        mkdir -p -v {wildcards.path}footprints/sm
        #
        mkdir -p -v {wildcards.path}footprints/ss/sites
        mkdir -p -v {wildcards.path}footprints/ss/insertions
        mkdir -p -v {wildcards.path}footprints/ss/stats
        mkdir -p -v {wildcards.path}footprints/ss/aggregated
        #
        mkdir -p -v {wildcards.path}footprints/sm/sites
        mkdir -p -v {wildcards.path}footprints/sm/insertions
        mkdir -p -v {wildcards.path}footprints/sm/stats
        mkdir -p -v {wildcards.path}footprints/sm/aggregated
        #
        touch {output}
        """


########################################################################################################################################
#### PREPROCESSING #####################################################################################################################
########################################################################################################################################

## This rule determines what is run in the preprocessing spooling option
rule AGGREGATOR_preprocessing:
    input:
        "{path}operations/directories/all_dirs.built",
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.peaks",
        #"{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.metrics",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw"
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.preprocessing"
    shell:
        """
        touch {output}
        """

## Gunzip the fastq files
rule STEP1_gunzip:
    input:
        a="{path}operations/directories/all_dirs.built",
        b="{path}fq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R{read}.fastq.gz"        
    output:
        "{path}preprocessing/1/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R{read}.fastq"
    shell:
        "gunzip -c {input.b} > {output}"

## Fastq QC filtering
rule STEP2_fastp_filtering:
    input:
        a="{path}preprocessing/1/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.fastq",
        b="{path}preprocessing/1/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.fastq"
    output:
        a="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        b="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    threads:
        6
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    conda:
    	"resources/envs/fastp.yaml"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.a} -O {output.b} -w {threads} -h {wildcards.path}metrics/fastq/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.quality.html -j {wildcards.path}metrics/fastq/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.quality.json"
  
## Align reads to reference genome
rule STEP3_refgenome_align:
    input:
        a="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        b="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4/{sample}.rep{repnum}.ref{refgenome}_L{lane}.sam"
    threads:
        12
    conda:
        "resources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "bowtie2 -q -p {threads} -X2000 -x genomes/{wildcards.refgenome}/{wildcards.refgenome} -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/align/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.al"

## Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP4_coordinate_sort_sam:
    ## A common error can occur here because samtools doesn't overwrite files,
    ## So if a run stopped while tmp output files were being sorted,
    ## Trying to sort again will cause an exception.
    ##  Delete the target files first to prevent this
    input:
        "{path}preprocessing/4/{sample}.rep{repnum}.ref{refgenome}_L{lane}.sam",
    output:
        "{path}preprocessing/5/{sample}.rep{repnum}.ref{refgenome}_L{lane}.cs.sam"
    conda:
        "resources/envs/samtools.yaml"
    shell:
        """
        rm -rf {wildcards.path}preprocessing/5/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}*
        samtools sort {input} -o {output} -O sam
        """

## Remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
rule STEP5_blacklist_filter_and_bam_conversion:
    input:
        "{path}preprocessing/5/{sample}.rep{repnum}.ref{refgenome}_L{lane}.cs.sam"
    output:
        a="{path}preprocessing/6/{sample}.rep{repnum}.ref{refgenome}_L{lane}.blacklist",
        b="{path}preprocessing/7/{sample}.rep{repnum}.ref{refgenome}_L{lane}.blrm.bam"
    conda:
        "resources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"
    
## Remove reads mapping to mitochondrial DNA
rule STEP6_chrM_contamination:
    input:
        "{path}preprocessing/7/{sample}.rep{repnum}.ref{refgenome}_L{lane}.blrm.bam"
    output:
        a="{path}preprocessing/8/{sample}.rep{repnum}.ref{refgenome}_L{lane}.mitochondrial",
        b="{path}preprocessing/9/{sample}.rep{repnum}.ref{refgenome}_L{lane}.good.bam"
    conda:
        "resources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/mtdna/mtdna.extents.bed -U {output.b} -@ 4 {input}"

## Add @RG tags to the reads and perform coordinate sorting
rule STEP7_add_rg_and_sort_bam:
    input:
        "{path}preprocessing/9/{sample}.rep{repnum}.ref{refgenome}_L{lane}.good.bam"
    output:
        "{path}preprocessing/a/{sample}.rep{repnum}.ref{refgenome}_L{lane}.tagged.bam"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 1
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.lane} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.lane}.{wildcards.sample} \
        RGSM={wildcards.sample}"
    
## Clean the bam file
rule STEP8_clean_bam:
    input:
        "{path}preprocessing/a/{sample}.rep{repnum}.ref{refgenome}_L{lane}.tagged.bam"
    output:
        "{path}preprocessing/b/{sample}.rep{repnum}.ref{refgenome}_L{lane}.clean.bam"
    conda:
        "resources/envs/picard.yaml"
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"
    
## Merge reads from different lanes
rule STEP9_merge_lanes:
    input:
        a="{path}preprocessing/b/{sample}.rep{repnum}.ref{refgenome}_L1.clean.bam",
        b="{path}preprocessing/b/{sample}.rep{repnum}.ref{refgenome}_L2.clean.bam",
        c="{path}preprocessing/b/{sample}.rep{repnum}.ref{refgenome}_L3.clean.bam",
        d="{path}preprocessing/b/{sample}.rep{repnum}.ref{refgenome}_L4.clean.bam"
    output:
        "{path}preprocessing/c/{sample}.rep{repnum}.ref{refgenome}.merged.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 1
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        I={input.d} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

## Remove PCR duplicate reads
rule STEP10_remove_pcr_duplicates:
    input:
        "{path}preprocessing/c/{sample}.rep{repnum}.ref{refgenome}.merged.bam"
    output:
        "{path}preprocessing/d/{sample}.rep{repnum}.ref{refgenome}.dp.bam"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 3
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output} \
        M={wildcards.path}metrics/duplication/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}.duplication \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
    
## Filter reads for only uniquely mapping
rule STEP11_mapq_filter:
    input:
        "{path}preprocessing/d/{sample}.rep{repnum}.ref{refgenome}.dp.bam"
    output:
        "{path}preprocessing/e/{sample}.rep{repnum}.ref{refgenome}.final.bam"
    conda:
        "resources/envs/samtools.yaml"
    resources:
        run_time=lambda params, attempt: attempt * 4
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
    
## Move bam to parent directory and rename
rule STEP12_move_bam:
    input:
        "{path}preprocessing/e/{sample}.rep{repnum}.ref{refgenome}.final.bam"
    output:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam"
    shell:
        """
        cp {wildcards.path}preprocessing/e/*rep{wildcards.repnum}*.bam {wildcards.path}aligned/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}.bam
        touch {output}
        """

## Build the .bai index for the processed bam file
rule STEP13_build_bai_index:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam"
    output:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 3
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"
    
## Make a bigwig file from the bam file
rule STEP14_make_bigwig:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw"
    conda:
        "resources/envs/deeptools.yaml"
    threads:
        10
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p {threads} -v"


########################################################################################################################################
#### PEAK CALLING ######################################################################################################################
########################################################################################################################################

## This rule determines what is done for peak calling
rule AGGREGATOR_peak_calling:
    input:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p01_peaks.narrowPeak",
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak",
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak",
        "{path}peaks/ln/{sample}.rep{repnum}.ref{refgenome}_localnorm_p01_peaks.narrowPeak",
        "{path}peaks/ln/{sample}.rep{repnum}.ref{refgenome}_localnorm_p001_peaks.narrowPeak",
        "{path}peaks/ln/{sample}.rep{repnum}.ref{refgenome}_localnorm_p0001_peaks.narrowPeak",
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak"
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.peaks"
    shell:
        "touch {output}"

## Call peaks using global normalization. pvalue 0.01
rule STEP16_MACS2_peaks_global_normilization_p01:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p01_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm_p01 --outdir {wildcards.path}peaks/gn --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

## Call peaks using global normalization. pvalue 0.001
rule STEP16_MACS2_peaks_global_normilization_p001:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm_p001 --outdir {wildcards.path}peaks/gn --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"

## Call peaks using global normalization. pvalue 0.0001
rule STEP16_MACS2_peaks_global_normilization_p0001:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm_p0001 --outdir {wildcards.path}peaks/gn --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.0001"

## Call peaks using local normalization. pvalue 0.01
rule STEP17_MACS2_peaks_local_normalization_p01:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/ln/{sample}.rep{repnum}.ref{refgenome}_localnorm_p01_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_localnorm_p01 --outdir {wildcards.path}peaks/ln --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

## Call peaks using local normalization. pvalue 0.001
rule STEP17_MACS2_peaks_local_normalization_p001:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/ln/{sample}.rep{repnum}.ref{refgenome}_localnorm_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_localnorm_p001 --outdir {wildcards.path}peaks/ln --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.001"

## Call peaks using local normalization. pvalue 0.001
rule STEP17_MACS2_peaks_local_normalization_p0001:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/ln/{sample}.rep{repnum}.ref{refgenome}_localnorm_p0001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_localnorm_p0001 --outdir {wildcards.path}peaks/ln --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.0001"


########################################################################################################################################
#### METRICS ###########################################################################################################################
########################################################################################################################################

## This rule determines what is done for the preprocessing metrics
rule AGGREGATOR_metrics:
    input:
        "{path}metrics/genomecov/{sample}.rep{repnum}.ref{refgenome}.peaks.txt",
        "{path}metrics/size/{sample}.rep{repnum}.ref{refgenome}.RData",
        "{path}metrics/myco/{sample}.rep{repnum}.ref{refgenome}_L1.alignment",
        "{path}metrics/fragsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_4.svg",
        "{path}metrics/peakideogram/{sample}.rep{repnum}.ref{refgenome}.svg",
        "{path}metrics/peakanno/{sample}.rep{repnum}.ref{refgenome}.annotate.peaks.complete"
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.metrics"
    shell:
        "touch {output}"

## Check for mycoplasma contamination
rule METRICS_mycoplasma_align:
    input:
        a="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        b="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    output:
        a="{path}preprocessing/3/{sample}.rep{repnum}.ref{refgenome}_L{lane}.myco",
        b="{path}metrics/myco/{sample}.rep{repnum}.ref{refgenome}_L{lane}.alignment"
    threads:
        6
    conda:
        "resources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    shell:
        "bowtie2 -q -p {threads} -X2000 -x genomes/myco/myco -1 {input.a} -2 {input.b} -S {output.a} 2>{output.b}"

## Calculate percent genome coverage from peaks with global normalization
rule METRICS_percent_peak_genome_coverage_globalnorm:
    input:
        a="{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p01_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/genomecov/{sample}.rep{repnum}.ref{refgenome}.peaks.txt"
    conda:
        "resources/envs/bedops.yaml"
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"

## Count the total number of reads in the sample
rule METRICS_sample_total_reads:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}metrics/size/{sample}.rep{repnum}.ref{refgenome}.RData"
    conda:
        "resources/envs/countSampleReads.yaml"
    script:
        "resources/scripts/preprocessing/countTotalSampleReads.R"

#### Generate the fragment size distribution graph
rule METRICS_fragment_size_distribution:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "resources/functions/atacFunctions.R"
    output:
        "{path}metrics/fragsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_1.svg",
        "{path}metrics/fragsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_2.svg",
        "{path}metrics/fragsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_3.svg",
        "{path}metrics/fragsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_4.svg"
    conda:
        "resources/envs/fragSizeDistribution.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/preprocessing/generateFragSizeDistribution.R"

#### Generate peak ideogram
rule METRICS_generate_peak_idiogram:
    input:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "{path}metrics/peakideogram/{sample}.rep{repnum}.ref{refgenome}.svg"
    conda:
        "resources/envs/generatePeakIdeogram.yaml"
    script:
        "resources/scripts/preprocessing/generatePeakIdeogram.R"

#### Generate peak annotations
rule METRICS_generate_peak_annotations:
    input:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "{path}metrics/peakanno/{sample}.rep{repnum}.ref{refgenome}.annotate.peaks.complete"
    conda:
        "resources/envs/annotatePeaks.yaml"
    script:
        "resources/scripts/preprocessing/annotatePeaks.R"


########################################################################################################################################
#### FOOTPRINTING ######################################################################################################################
########################################################################################################################################

## Step 1, determine the binding sites
rule FOOTPRINTING_generate_binding_sites_sm:
    input:
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/generateBindingSites.R"

## Step 2, generate the insertion matrix
rule FOOTPRINTING_generate_insertion_matrix_sm:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}footprints/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.bind.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/sm/insertions/{sample}.rep{repnum}.ref{refgenome}.{gene}.complete"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateInsertionMatrix.R"

## Step 3, generate the footprint stats table
rule FOOTPRINTING_generate_footprint_stats_sm:
    input:
        "{path}footprints/sm/insertions/{sample}.rep{repnum}.ref{refgenome}.{gene}.complete",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.stats.complete"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateFootprintStatsByChr.R"

## Step 4, aggregate the by chr stats tables into a single file
rule FOOTPRINTING_aggregate_footprint_stats_sm:
    input:
        "{path}footprints/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.stats.complete",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/sm/aggregated/{sample}.rep{repnum}.ref{refgenome}.{gene}.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/aggregateFootprintStats.R"

##
rule AGGREGATOR_footprinting:
    input:
        expand("{{path}}footprints/sm/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.RData", genename=config["geneNames"])
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.uncorrected"
    shell:
        "touch {output}"

########################################################################################################################################
#### FOOTPRINTING - JASPAR #############################################################################################################
########################################################################################################################################

## Step 1, determine the binding sites
rule FOOTPRINTING_generate_binding_sites_sm_JASPAR:
    input:
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "resources/functions/atacFunctions.R",
        "resources/jaspar/pwm/{ID}.pwm.RData"        
    output:
        "{path}footprints/jas/sites/{sample}.rep{repnum}.ref{refgenome}.{ID}.bind.jp.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/generateBindingSitesJaspar.R"

## Step 2, generate the insertion matrix
rule FOOTPRINTING_generate_insertion_matrix_sm_JASPAR:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}footprints/jas/sites/{sample}.rep{repnum}.ref{refgenome}.{ID}.bind.jp.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/jas/insertions/{sample}.rep{repnum}.ref{refgenome}.{ID}.complete.jp"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateInsertionMatrixJaspar.R"

## Step 3, generate the footprint stats table
rule FOOTPRINTING_generate_footprint_stats_sm_JASPAR:
    input:
        "{path}footprints/jas/insertions/{sample}.rep{repnum}.ref{refgenome}.{ID}.complete.jp",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/jas/stats/{sample}.rep{repnum}.ref{refgenome}.{ID}.stats.complete.jp"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateFootprintStatsByChrJaspar.R"

## Step 4, aggregate the by chr stats tables into a single file
rule FOOTPRINTING_aggregate_footprint_stats_sm_JASPAR:
    input:
        "{path}footprints/jas/stats/{sample}.rep{repnum}.ref{refgenome}.{ID}.stats.complete.jp",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/jas/aggregated/{sample}.rep{repnum}.ref{refgenome}.{ID}.RData.jp"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/aggregateFootprintStatsJaspar.R"

##
rule AGGREGATOR_footprinting_JASPAR:
    input:
        expand("{{path}}footprints/jas/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{jasparID}.RData.jp", jasparID=config["jasparNames"])
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.jp"
    shell:
        "touch {output}"


########################################################################################################################################
#### FOOTPRINTING - FORKED #############################################################################################################
########################################################################################################################################

## Step 1, determine the binding sites
rule FOOTPRINTING_generate_binding_sites_sm_forked:
    input:
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/generateBindingSites.R"

## Step 2, generate the insertion matrix
rule FOOTPRINTING_generate_insertion_matrix_sm_forked:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}footprints/forked/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.bind.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.{chr}.ins.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateInsertionMatrixForked.R"

## Step 3, generate the footprint stats table
rule FOOTPRINTING_generate_footprint_stats_sm_forked:
    input:
        "{path}footprints/forked/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.{chr}.ins.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.{chr}.stats.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateFootprintStatsForked.R"

rule FOOTPRINTING_fork_footprint_stats_forked:
    input:
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr1.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr2.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr3.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr4.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr5.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr6.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr7.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr8.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr9.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr10.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr11.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr12.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr13.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr14.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr15.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr16.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr17.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr18.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr19.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr20.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr21.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr22.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chrX.stats.RData",
        "{path}footprints/forked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chrY.stats.RData"
    output:
        "{path}footprints/forked/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.complete"
    shell:
        "touch {output}"

## Step 4, aggregate the by chr stats tables into a single file
rule FOOTPRINTING_aggregate_footprint_stats_sm_forked:
    input:
        "{path}footprints/forked/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.complete",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/sm/aggregated/{sample}.rep{repnum}.ref{refgenome}.{gene}.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/aggregateFootprintStatsForked.R"

##
rule AGGREGATOR_footprinting_forked:
    input:
        expand("{{path}}footprints/forked/sm/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.RData", genename=config["geneNames"])
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.uncorrected.forked"
    shell:
        "touch {output}"

########################################################################################################################################
#### FOOTPRINTING - FORKED JASPAR ######################################################################################################
########################################################################################################################################

## Step 1, determine the binding sites
rule FOOTPRINTING_generate_binding_sites_sm_forked_jaspar:
    input:
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/jas/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.bind.jp.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/generateBindingSitesJaspar.R"

## Step 2, generate the insertion matrix
rule FOOTPRINTING_generate_insertion_matrix_sm_forked_jaspar:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}footprints/forked/jas/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.bind.jp.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/jas/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.{chr}.ins.jp.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateInsertionMatrixForked.R"

## Step 3, generate the footprint stats table
rule FOOTPRINTING_generate_footprint_stats_sm_forked_jaspar:
    input:
        "{path}footprints/forked/jas/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.{chr}.ins.jp.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.{chr}.stats.jp.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
    script:
        "resources/scripts/footprinting/generateFootprintStatsForked.R"

rule FOOTPRINTING_fork_footprint_stats_forked_jaspar:
    input:
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr1.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr2.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr3.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr4.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr5.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr6.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr7.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr8.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr9.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr10.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr11.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr12.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr13.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr14.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr15.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr16.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr17.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr18.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr19.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr20.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr21.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chr22.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chrX.stats.jp.RData",
        "{path}footprints/forked/jas/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.chrY.stats.jp.RData"
    output:
        "{path}footprints/forked/jas/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.jp.complete"
    shell:
        "touch {output}"

## Step 4, aggregate the by chr stats tables into a single file
rule FOOTPRINTING_aggregate_footprint_stats_sm_forked_jaspar:
    input:
        "{path}footprints/forked/jas/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.jp.complete",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/forked/jas/sm/aggregated/{sample}.rep{repnum}.ref{refgenome}.{gene}.jp.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/aggregateFootprintStatsForked.R"

##
rule AGGREGATOR_footprinting_forked_jaspar:
    input:
        expand("{{path}}footprints/forked/jas/sm/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.jp.RData", genename=config["geneNames"])
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.jp.forked"
    shell:
        "touch {output}"


########################################################################################################################################
#### FOOTPRINTING - CHUNKED ############################################################################################################
########################################################################################################################################

## Step 1, determine the binding sites
rule FOOTPRINTING_generate_binding_sites_sm_chunked:
    input:
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/chunked/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/generateBindingSites.R"

## Step 2, generate the insertion matrix
rule FOOTPRINTING_generate_insertion_matrix_sm_chunked:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}footprints/chunked/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.cnk.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/chunked/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.ins.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 6
    script:
        "resources/scripts/footprinting/generateInsertionMatrixChunked.R"

## Step 3, generate the footprint stats table
rule FOOTPRINTING_generate_footprint_stats_sm_chunked:
    input:
        "{path}footprints/chunked/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.ins.cnk.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.stats.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 6
    script:
        "resources/scripts/footprinting/generateFootprintStatsChunked.R"

rule FOOTPRINTING_fork_footprint_stats_chunked:
    input:
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.1.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.2.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.3.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.4.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.5.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.6.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.7.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.8.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.9.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.10.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.11.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.12.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.13.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.14.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.15.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.16.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.17.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.18.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.19.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.20.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.21.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.22.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.23.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.24.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.25.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.26.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.27.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.28.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.29.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.30.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.31.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.32.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.33.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.34.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.35.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.36.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.37.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.38.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.39.stats.cnk.RData",
        "{path}footprints/chunked/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.40.stats.cnk.RData"
    output:
        "{path}footprints/chunked/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.cnk.complete"
    shell:
        "touch {output}"

## Step 4, aggregate the by chr stats tables into a single file
rule FOOTPRINTING_aggregate_footprint_stats_sm_chunked:
    input:
        "{path}footprints/chunked/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.cnk.complete",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/chunked/sm/aggregated/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/aggregateFootprintStatsChunked.R"

##
rule AGGREGATOR_footprinting_chunked:
    input:
        expand("{{path}}footprints/chunked/sm/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.cnk.RData", genename=config["geneNames"])
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.pwm{matchScore}.uncorrected.chunked"
    shell:
        "touch {output}"


########################################################################################################################################
#### FOOTPRINTING - JASPAR CHUNKED #####################################################################################################
########################################################################################################################################

## Step 1, determine the binding sites
rule FOOTPRINTING_generate_binding_sites_sm_JASPAR_chunked:
    input:
        "{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "resources/functions/atacFunctions.R",
        "resources/jaspar/pwm/{ID}.pwm.RData"        
    output:
        "{path}footprints/jaspar/chunked/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.bind.jp.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/generateBindingSitesJasparChunked.R"

## Step 2, generate the insertion matrix
rule FOOTPRINTING_generate_insertion_matrix_sm_jaspar_chunked:
    input:
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}footprints/jaspar/chunked/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.bind.jp.cnk.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/jaspar/chunked/sm/insertions/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.{chunk}.ins.jp.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateInsertionMatrixChunkedJaspar.R"

## Step 3, generate the footprint stats table
rule FOOTPRINTING_generate_footprint_stats_sm_jaspar_chunked:
    input:
        "{path}footprints/jaspar/chunked/sm/insertions/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.{chunk}.ins.jp.cnk.RData",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.{chunk}.stats.jp.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateFootprintStatsChunkedJaspar.R"

rule FOOTPRINTING_fork_footprint_stats_jaspar_chunked:
    input:
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.1.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.2.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.3.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.4.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.5.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.6.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.7.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.8.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.9.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.10.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.11.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.12.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.13.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.14.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.15.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.16.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.17.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.18.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.19.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.20.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.21.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.22.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.23.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.24.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.25.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.26.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.27.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.28.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.29.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.30.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.31.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.32.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.33.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.34.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.35.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.36.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.37.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.38.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.39.stats.jp.cnk.RData",
        "{path}footprints/jaspar/chunked/sm/stats/{ID}/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.40.stats.jp.cnk.RData"
    output:
        "{path}footprints/jaspar/chunked/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.jp.cnk.complete"
    shell:
        "touch {output}"

## Step 4, aggregate the by chr stats tables into a single file
rule FOOTPRINTING_aggregate_footprint_stats_sm_jaspar_chunked:
    input:
        "{path}footprints/jaspar/chunked/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.jp.cnk.complete",
        "resources/functions/atacFunctions.R"
    output:
        "{path}footprints/jaspar/chunked/sm/aggregated/{sample}.rep{repnum}.ref{refgenome}.{ID}.pwm{matchScore}.jp.cnk.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "resources/scripts/footprinting/aggregateFootprintStatsChunkedJaspar.R"

##
rule AGGREGATOR_footprinting_jaspar_chunked:
    input:
        expand("{{path}}footprints/jaspar/chunked/sm/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{ID}.pwm{{matchScore}}.jp.cnk.RData", ID=config["jasparNames"])
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.pwm{matchScore}.jp.uncorrected.chunked"
    shell:
        "touch {output}"



########################################################################################################################################
#### SEQBIAS MODELS ####################################################################################################################
########################################################################################################################################

## This rule determines what is run for the footprinting analysis
rule AGGREGATOR_seqbias_models:
    input:
        "{path}seqbias/{sample}.rep{repnum}.ref{refgenome}.seqbias_sample_reads_against_sample_peaks.yml",
        "{path}seqbias/{sample}.rep{repnum}.ref{refgenome}.seqbias_sample_reads_against_hg38_chr1.yml"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.seqbias_models.complete"
    shell:
        "touch {output}"

##
rule SEQBIAS_write_sample_peaks_to_BED:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "snakeResources/functions/atacFunctions.R"
    output:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.bed"
    conda:
        "snakeResources/envs/modelSeqbias.yaml"
    script:
        "snakeResources/scripts/convertPeaksToBEDinterval.R"


########################################################################################################################################
#### SAMPLE SPECIFIC MERGE PEAKS #######################################################################################################
########################################################################################################################################

rule MERGE_sample_peaks_lncap:
    input:
        "data/pros/lncap/cr01/peaks/gn/LNCaP-CR-01.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/cr02/peaks/gn/LNCaP-CR-02.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/cr04/peaks/gn/LNCaP-CR-04.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/cr05/peaks/gn/LNCaP-CR-05.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/cr07/peaks/gn/LNCaP-CR-07.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/cr08/peaks/gn/LNCaP-CR-08.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/wt01/peaks/gn/LNCaP-WT-01.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/pros/lncap/wt02/peaks/gn/LNCaP-WT-02.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "data/pros/lncap/cr01/peaks/sm/LNCaP-CR-01.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr02/peaks/sm/LNCaP-CR-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr04/peaks/sm/LNCaP-CR-04.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr05/peaks/sm/LNCaP-CR-05.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr07/peaks/sm/LNCaP-CR-07.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr08/peaks/sm/LNCaP-CR-08.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/wt01/peaks/sm/LNCaP-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/wt02/peaks/sm/LNCaP-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "resources/envs/mergeSamplePeaks.yaml"
    script:
        "resources/scripts/preprocessing/mergeSamplePeaks.R"

rule MERGE_sample_peaks_cosma:
    input:
        "data/cosma/DAbaz2b/peaks/gn/DonorA_Baz2B.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/cosma/DAluf/peaks/gn/DonorA_Luf.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/cosma/DAprog/peaks/gn/DonorA_Progenitor.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/cosma/DBbaz2b/peaks/gn/DonorB_Baz2B.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/cosma/DBluf/peaks/gn/DonorB_Luf.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/cosma/DBprog/peaks/gn/DonorB_Progenitor.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "data/cosma/UNd/peaks/gn/Undetermined.rep1.refhg38_globalnorm_p01_peaks.narrowPeak",
        "resources/functions/atacFunctions.R"
    output:
        "data/cosma/DAbaz2b/peaks/sm/DonorA_Baz2B.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/cosma/DAluf/peaks/sm/DonorA_Luf.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/cosma/DAprog/peaks/sm/DonorA_Progenitor.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/cosma/DBbaz2b/peaks/sm/DonorB_Baz2B.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/cosma/DBluf/peaks/sm/DonorB_Luf.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/cosma/DBprog/peaks/sm/DonorB_Progenitor.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/cosma/UNd/peaks/sm/Undetermined.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "resources/envs/mergeSamplePeaks.yaml"
    script:
        "resources/scripts/preprocessing/mergeSamplePeaks.R"
