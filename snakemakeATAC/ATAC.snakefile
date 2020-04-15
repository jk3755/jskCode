########################################################################################################################################
#### IMPORT CONFIGURATION FILE #########################################################################################################
########################################################################################################################################
configfile: "resources/config/config.yaml"


########################################################################################################################################
#### IMPORT SPOOLING MODULES ###########################################################################################################
########################################################################################################################################
include: "resources/modules/spool_LNCAP.snakefile"
include: "resources/modules/spool_CICCIA.snakefile"
include: "resources/modules/spool_COSMA.snakefile"


########################################################################################################################################
#### FULL ANALYSIS #####################################################################################################################
########################################################################################################################################
## This rule determines what is run in the full analysis spooling option
## Also cleans up the intermediate preprocessing data
rule AGGREGATOR_full_analysis:
    input:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.preprocessing"
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
## Generate the required directories, some software may fail when trying to create files if directory doesnt exist
## Important to have separated directory structure to avoid dirs with too many files and confusing snakemake metadata
## Avoid repeating words in the full paths to intermediate files

# Directory generation aggregator rule
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

# The parent directory
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

# Flag files for marking when a step is complete in the pipeline
rule DIR_operations:
    output:
        "{path}operations/directories/operations_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations/directories
        mkdir -p -v {wildcards.path}operations/aggregators
        touch {output}
        """

# Directories for preprocessing data
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

# For QC metrics files
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

# For peak calling
rule DIR_peaks:
    output:
        "{path}operations/directories/peaks_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}peaks/ln {wildcards.path}peaks/gn {wildcards.path}peaks/sm
        touch {output}
        """

# For footprinting
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
        #"{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.peaks",
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
 
# Using bowtie2, may deprecate because STAR gets higher alignment rates
# Align reads to reference genome
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

# ## Using STAR
# ## Align reads to reference genome
# rule STEP3_refgenome_align:
#     input:
#         a="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
#         b="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
#     output:
#         "{path}preprocessing/4/{sample}.rep{repnum}.ref{refgenome}_L{lane}.sam"
#     threads:
#         15
#     conda:
#         "resources/envs/star.yaml"
#     resources:
#         mem_mb=lambda params, attempt: attempt * 25000,
#         run_time=lambda params, attempt: attempt * 12
#     shell:
#         "STAR --genomeDir genomes/star/hg38 \
# 		--outFileNamePrefix test/test \
# 		--readFilesIn Merged_DonorA_BAZ2B_R1.fastq.gz Merged_DonorA_BAZ2B_R2.fastq.gz \
# 		--readFilesCommand zcat \
# 		--runThreadN 20 \
# 		--outFilterMultimapScoreRange 1 \
# 		--outFilterMultimapNmax 1 \
# 		--outFilterMismatchNoverReadLmax 0.1 \
# 		--alignIntronMax 1 \
# 		--alignMatesGapMax 1000 \
# 		--alignEndsType EndToEnd \
# 		--genomeLoad NoSharedMemory \
# 		--seedSearchStartLmax 30 \
# 		--limitBAMsortRAM 0 \
# 		--outSAMattributes All \
# 		--outSAMunmapped None \
# 		--outSAMtype BAM SortedByCoordinate \
# 		--outSAMheaderHD @HD VN:1.4"
		
## Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP4_coordinate_sort_sam:
    ## A common error can occur here because samtools doesn't overwrite files,
    ## So if a run stopped while tmp output files were being sorted,
    ## Trying to sort again will cause an exception.
    ## Delete the target files first to prevent this
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
        ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam")
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
        a=ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam"),
        b=ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai")
    output:
        "{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw"
    conda:
        "resources/envs/deeptools.yaml"
    benchmark:
        "{path}benchmark/bigwig/{sample}.rep{repnum}.ref{refgenome}.bigwig.txt"
    threads:
        10
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
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
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
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
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm --outdir {wildcards.path}peaks/gn --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.0001"

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
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
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
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
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
rule FOOTPRINTING_generate_binding_sites_sm:
    input:
        ancient("{path}peaks/sm/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateBindingSites.R"

rule FOOTPRINTING_generate_insertion_matrix_sm:
    input:
        ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam"),
        ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"),
        ancient("{path}footprints/sm/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.ins.totalchunk{totalchunk}.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateInsertionMatrix.R"

rule FOOTPRINTING_generate_footprint_stats_sm:
    input:
        ancient("{path}footprints/sm/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.ins.totalchunk{totalchunk}.RData"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/sm/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.stats.totalchunk{totalchunk}.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateFootprintStats.R"

rule FOOTPRINTING_fork_footprint_stats_sm:
    input:
        ancient(expand("{{path}}footprints/sm/stats/{{gene}}/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{{gene}}.pwm{{matchScore}}.{chunkNum}.stats.totalchunk{{totalchunk}}.RData", chunkNum=config["chunkFootprinting"]))
    output:
        "{path}footprints/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.stats.totalchunk{totalchunk}.complete"
    shell:
        "touch {output}"

rule FOOTPRINTING_aggregate_footprint_stats_sm:
    input:
        ancient("{path}footprints/sm/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.stats.totalchunk{totalchunk}.complete"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/sm/aggregated/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.sm.totalchunk{totalchunk}.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/aggregateFootprintStats.R"

rule AGGREGATOR_footprinting:
    input:
        ancient(expand("{{path}}footprints/sm/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.sm.totalchunk{{totalchunk}}.RData", genename=config["nfyGenes"]))
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.pwm{matchScore}.uncorrected.sm.totalchunk{totalchunk}"
    shell:
        "touch {output}"

########################################################################################################################################
#### FOOTPRINTING - NO PEAKS ###########################################################################################################
########################################################################################################################################
rule FOOTPRINTING_generate_binding_sites_nopeak:
    input:
        ancient("/ifs/scratch/c2b2/ac_lab/jk3755/atac/resources/hg38.bed"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/nopeak/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.nopeak.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateBindingSites.R"

rule FOOTPRINTING_generate_insertion_matrix_nopeak:
    input:
        ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam"),
        ancient("{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"),
        ancient("{path}footprints/nopeak/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.nopeak.RData"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/nopeak/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.ins.nopeak.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateInsertionMatrix.R"

rule FOOTPRINTING_generate_footprint_stats_nopeak:
    input:
        ancient("{path}footprints/nopeak/insertions/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.ins.nopeak.RData"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/nopeak/stats/{gene}/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.{chunk}.stats.nopeak.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/generateFootprintStats.R"

rule FOOTPRINTING_fork_footprint_stats_nopeak:
    input:
        ancient(expand("{{path}}footprints/nopeak/stats/{{gene}}/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{{gene}}.pwm{{matchScore}}.{chunkNum}.stats.nopeak.RData", chunkNum=config["chunkFootprinting"]))
    output:
        "{path}footprints/nopeak/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.stats.nopeak.complete"
    shell:
        "touch {output}"

rule FOOTPRINTING_aggregate_footprint_stats_nopeak:
    input:
        ancient("{path}footprints/nopeak/stats/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.stats.nopeak.complete"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "{path}footprints/nopeak/aggregated/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.nopeak.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    script:
        "resources/scripts/footprinting/aggregateFootprintStats.R"

rule AGGREGATOR_footprinting_nopeak:
    input:
        ancient(expand("{{path}}footprints/nopeak/aggregated/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.nopeak.RData", genename=config["nfyGenes"]))
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.pwm{matchScore}.nopeak"
    shell:
        "touch {output}"


########################################################################################################################################
#### PEAKS - COSMA #####################################################################################################################
########################################################################################################################################

## pvalue 0.05 ##
rule cosma_MACS2_peaks_global_normilization_p005_1:
    input:
        a="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_merged_sorted.bam",
        b="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_merged_sorted.bam.bai"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_global_p005_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n baz2b_global_p005 --outdir /ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.005"

rule cosma_MACS2_peaks_global_normilization_p005_2:
    input:
        a="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_merged_sorted.bam",
        b="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_merged_sorted.bam.bai"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_global_p005_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n luciferase_global_p005 --outdir /ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.005"

rule cosma_MACS2_peaks_global_normilization_p005_3:
    input:
        a="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_merged_sorted.bam",
        b="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_merged_sorted.bam.bai"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_global_p005_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n progenitor_global_p005 --outdir /ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.005"

rule target_cosma_peaks_p005:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_global_p005_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_global_p005_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_global_p005_peaks.narrowPeak"

## pvalue 0.01 ##
rule cosma_MACS2_peaks_global_normilization_p001_1:
    input:
        a="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_merged_sorted.bam",
        b="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_merged_sorted.bam.bai"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_global_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n baz2b_global_p001 --outdir /ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"


rule cosma_MACS2_peaks_global_normilization_p001_2:
    input:
        a="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_merged_sorted.bam",
        b="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_merged_sorted.bam.bai"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_global_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n luciferase_global_p001 --outdir /ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"


rule cosma_MACS2_peaks_global_normilization_p001_3:
    input:
        a="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_merged_sorted.bam",
        b="/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_merged_sorted.bam.bai"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_global_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n progenitor_global_p001 --outdir /ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"

rule target_cosma_peaks_p001:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/progenitor_global_p001_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/luciferase_global_p001_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/cosma_merged/baz2b_global_p001_peaks.narrowPeak"

########################################################################################################################################
#### Count binding sites for all genes at 80% ciccia ###################################################################################
########################################################################################################################################

rule count_sites_ciccia_80:
    input:
        ancient("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/CTRL_UNTR_R1/peaks/sm/MDA_MB436_CTRL_UNTR.rep1.refhg38_sample_merged_peaks.narrowPeak"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/80/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 10
    script:
        "resources/scripts/footprinting/countBindingSites.R"

rule AGGREGATOR_count_sites_ciccia_80:
    input:
        ancient(expand("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/80/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.bind.RData", genename=config["allGenes"]))
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/80/{sample}.rep{repnum}.ref{refgenome}.pwm{matchScore}_finished.txt"
    shell:
        "touch {output}"
##################################################################################################################################################################################
rule count_sites_ciccia_85:
    input:
        ancient("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/CTRL_UNTR_R1/peaks/sm/MDA_MB436_CTRL_UNTR.rep1.refhg38_sample_merged_peaks.narrowPeak"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/85/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 10
    script:
        "resources/scripts/footprinting/countBindingSites.R"

rule AGGREGATOR_count_sites_ciccia_85:
    input:
        ancient(expand("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/85/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.bind.RData", genename=config["allGenes"]))
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/85/{sample}.rep{repnum}.ref{refgenome}.pwm{matchScore}_finished.txt"
    shell:
        "touch {output}"

##################################################################################################################################################################################
rule count_sites_ciccia_90:
    input:
        ancient("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/CTRL_UNTR_R1/peaks/sm/MDA_MB436_CTRL_UNTR.rep1.refhg38_sample_merged_peaks.narrowPeak"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/90/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 10
    script:
        "resources/scripts/footprinting/countBindingSites.R"

rule AGGREGATOR_count_sites_ciccia_90:
    input:
        ancient(expand("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/90/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.bind.RData", genename=config["allGenes"]))
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/90/{sample}.rep{repnum}.ref{refgenome}.pwm{matchScore}_finished.txt"
    shell:
        "touch {output}"

##################################################################################################################################################################################
rule count_sites_ciccia_95:
    input:
        ancient("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/CTRL_UNTR_R1/peaks/sm/MDA_MB436_CTRL_UNTR.rep1.refhg38_sample_merged_peaks.narrowPeak"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/95/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 10
    script:
        "resources/scripts/footprinting/countBindingSites.R"

rule AGGREGATOR_count_sites_ciccia_95:
    input:
        ancient(expand("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/95/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.bind.RData", genename=config["allGenes"]))
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/95/{sample}.rep{repnum}.ref{refgenome}.pwm{matchScore}_finished.txt"
    shell:
        "touch {output}"
##################################################################################################################################################################################
rule count_sites_ciccia_99:
    input:
        ancient("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/CTRL_UNTR_R1/peaks/sm/MDA_MB436_CTRL_UNTR.rep1.refhg38_sample_merged_peaks.narrowPeak"),
        ancient("resources/functions/atacFunctions.R")
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/99/{sample}.rep{repnum}.ref{refgenome}.{gene}.pwm{matchScore}.bind.RData"
    conda:
        "resources/envs/footprintAnalysis.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 10
    script:
        "resources/scripts/footprinting/countBindingSites.R"

rule AGGREGATOR_count_sites_ciccia_99:
    input:
        ancient(expand("/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/99/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.pwm{{matchScore}}.bind.RData", genename=config["allGenes"]))
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/99/{sample}.rep{repnum}.ref{refgenome}.pwm{matchScore}_finished.txt"
    shell:
        "touch {output}"

####
rule count_binding_sites_ciccia_all:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/80/ciccia.rep1.refhg38.pwm80_finished.txt",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/85/ciccia.rep1.refhg38.pwm85_finished.txt",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/90/ciccia.rep1.refhg38.pwm90_finished.txt",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/95/ciccia.rep1.refhg38.pwm95_finished.txt",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/99/ciccia.rep1.refhg38.pwm99_finished.txt"

##
rule count_binding_sites_ciccia_80:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/80/ciccia.rep1.refhg38.pwm80_finished.txt"

rule count_binding_sites_ciccia_85:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/85/ciccia.rep1.refhg38.pwm85_finished.txt"

rule count_binding_sites_ciccia_90:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/90/ciccia.rep1.refhg38.pwm90_finished.txt"

rule count_binding_sites_ciccia_95:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/95/ciccia.rep1.refhg38.pwm95_finished.txt"

rule count_binding_sites_ciccia_99:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/ciccia/count_sites/99/ciccia.rep1.refhg38.pwm99_finished.txt"

########################################################################################################################################
#### LNCAP - CPM NORMALIZE BAM #########################################################################################################
########################################################################################################################################

##
rule CPM_NORMALIZE_BAM:
    input:
        ancient("{path}aligned/{sample}.bam")
    output:
        "{path}cpm_normalized/{sample}.cpm_normalized.bw"
    conda:
        "resources/envs/deeptools.yaml"
    threads:
        4
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM -b {input} -o MDA_MB436_CTRL_UNTR_rep1_hg38_sort.bam.bw"

rule CPM_normalize_lncap:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/cr01/cpm_normalized/LNCaP-CR-01.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/cr02/cpm_normalized/LNCaP-CR-02.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/cr04/cpm_normalized/LNCaP-CR-04.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/cr05/cpm_normalized/LNCaP-CR-05.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/cr07/cpm_normalized/LNCaP-CR-07.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/cr08/cpm_normalized/LNCaP-CR-08.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/wt01/cpm_normalized/LNCaP-WT-01.rep1.refhg38.cpm_normalized.bw",
        "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/lncap/wt02/cpm_normalized/LNCaP-WT-02.rep1.refhg38.cpm_normalized.bw"

########################################################################################################################################
#### TCGA BREAST CANCER PEAK CALLING ###################################################################################################
########################################################################################################################################

##############################################################################
rule TCGA_MACS2_PEAK_CALLING:
    input:
        "/ifs/scratch/c2b2/ac_lab/av2729/BRCA-TCGA-ATAC-seq/BRCA_bedGraph/{sample}.insertions.bw.bedGraph"
    output:
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/{sample}.globalnorm_c0.0_l200_g30_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    shell:
        "macs2 bdgpeakcall -i {input} --o-prefix {wildcards.sample}.globalnorm --outdir /ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks"

####
rule TCGA_peak_calls_all:
    input:
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_000CFD9F_ADDF_4304_9E60_6041549E189C_X017_S06_L011_B1_T1_P040.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_000CFD9F_ADDF_4304_9E60_6041549E189C_X017_S06_L012_B1_T2_P046.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1B6783CB_7D24_4D13_908A_19CCAD4CFF34_X002_S06_L011_B1_T1_P001.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1B6783CB_7D24_4D13_908A_19CCAD4CFF34_X002_S06_L012_B1_T2_P001.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1BD15E25_9CE5_49AB_ADCF_41F857F3468D_X021_S04_L030_B1_T1_P054.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1BD15E25_9CE5_49AB_ADCF_41F857F3468D_X021_S04_L031_B1_T2_P053.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1D939DC3_EF0C_40BF_BC60_8C5D46345265_X021_S01_L024_B1_T1_P051.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1D939DC3_EF0C_40BF_BC60_8C5D46345265_X021_S01_L025_B1_T2_P052.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1FF8A889_24D2_4439_8BB4_C9AB39B7EE23_X005_S12_L028_B1_T1_P011.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1FF8A889_24D2_4439_8BB4_C9AB39B7EE23_X005_S12_L029_B1_T2_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2A65DC63_F8CC_4EF4_AB23_3F5FD880FB5E_X018_S08_L038_B1_T1_P044.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2A65DC63_F8CC_4EF4_AB23_3F5FD880FB5E_X018_S08_L039_B1_T2_P043.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2C46607F_72A2_4B40_A311_E1F29C59DCFA_X010_S01_L049_B1_T1_P016.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2C46607F_72A2_4B40_A311_E1F29C59DCFA_X010_S01_L050_B1_T2_P019.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2D277FFC_F83F_4882_98D7_525388B79C3B_X023_S01_L073_B1_T1_P057.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2D277FFC_F83F_4882_98D7_525388B79C3B_X023_S01_L074_B1_T2_P056.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_4DE199D4_0892_462C_BE6E_09413F7EF76E_X018_S03_L029_B1_T1_P045.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_4F06A0BF_EB31_45EE_BB09_5601502D86FE_X021_S07_L036_B1_T1_P052.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_4F06A0BF_EB31_45EE_BB09_5601502D86FE_X021_S07_L037_B1_T2_P051.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_5C54B79C_DA02_4B22_9FC2_3D61BFFC5559_X019_S11_L068_B1_T1_P046.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_5C54B79C_DA02_4B22_9FC2_3D61BFFC5559_X019_S11_L069_B1_T2_P046.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_5C60C974_536D_4BD9_ADC5_606CA5E60C04_X020_S10_L019_B1_T1_P047.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_6D0DAA1B_2ACC_49CC_90C5_07FC4260432E_X011_S05_L009_B1_T1_P026.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_6D0DAA1B_2ACC_49CC_90C5_07FC4260432E_X011_S05_L010_B1_T2_P025.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_6F22B7DA_85CA_4E9C_93A3_859878775DDB_X005_S06_L014_B1_T3_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_6F22B7DA_85CA_4E9C_93A3_859878775DDB_X005_S06_L015_B1_T4_P014.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_7C6A3AE4_E2EA_42B3_B3F1_81C19E6F2170_X023_S05_L081_B1_T1_P057.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_7C6A3AE4_E2EA_42B3_B3F1_81C19E6F2170_X023_S05_L082_B1_T2_P055.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_7D3BBF34_47F6_41EE_B5BF_6571C18984A8_X020_S03_L005_B1_T1_P053.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_7D3BBF34_47F6_41EE_B5BF_6571C18984A8_X020_S03_L006_B1_T2_P054.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8D1E6006_85CB_484A_8B5C_30766D90137B_X012_S03_L029_B1_T1_P024.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8D1E6006_85CB_484A_8B5C_30766D90137B_X012_S03_L030_B1_T2_P026.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8D51DE4E_9B6A_4C2D_BD11_5F66C2173002_X015_S06_L035_B1_T1_P035.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8D51DE4E_9B6A_4C2D_BD11_5F66C2173002_X015_S06_L036_B1_T2_P036.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8D66AC5D_CA8B_47FC_8C61_376A75AEB292_X019_S05_L056_B1_T1_P043.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8D66AC5D_CA8B_47FC_8C61_376A75AEB292_X019_S05_L057_B1_T2_P038.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8DB1A244_D21F_4A52_A70F_317C0CF8B4B7_X014_S11_L021_B1_T1_P033.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8DB1A244_D21F_4A52_A70F_317C0CF8B4B7_X014_S11_L022_B1_T2_P034.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_11AAB1D2_D77B_42AA_AA32_786A597B9797_X003_S06_L029_B1_T1_P002.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_11AAB1D2_D77B_42AA_AA32_786A597B9797_X003_S06_L030_B1_T2_P002.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_11D015DD_1250_48BC_8B5D_3262C97F164B_X012_S08_L039_B1_T1_P029.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_11D015DD_1250_48BC_8B5D_3262C97F164B_X012_S08_L040_B1_T2_P022.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_14AD76EE_12F9_40B3_8DCD_4A256E02CF8D_X020_S12_L022_B1_T1_P053.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_14AD76EE_12F9_40B3_8DCD_4A256E02CF8D_X020_S12_L023_B1_T2_P054.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_36EA9779_3126_45DE_8C5D_3E6B3680086C_X012_S06_L035_B1_T1_P022.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_36EA9779_3126_45DE_8C5D_3E6B3680086C_X012_S06_L036_B1_T2_P025.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_64BD3CE0_7385_461B_9870_AD6611DD886A_X014_S01_L001_B1_T1_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_64BD3CE0_7385_461B_9870_AD6611DD886A_X014_S01_L002_B1_T2_P037.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_74CCA3AC_8984_4207_AD8C_979E7596A5DC_X005_S02_L003_B1_T1_P014.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_74CCA3AC_8984_4207_AD8C_979E7596A5DC_X005_S02_L004_B1_T2_P015.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_78F79CE9_210D_47BF_B6AB_AA5327778FBE_X010_S03_L053_B1_T1_P021.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_78F79CE9_210D_47BF_B6AB_AA5327778FBE_X010_S03_L054_B1_T2_P019.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_85C6EA5B_FA06_416B_9D98_80339B9B03DC_X005_S04_L007_B1_T1_P015.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_86DA5538_752C_4F30_A411_0970F65BD115_X012_S07_L037_B1_T1_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_86DA5538_752C_4F30_A411_0970F65BD115_X012_S07_L038_B1_T2_P030.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_94AF19F0_1F2A_41EC_8CB6_96C76227811F_X013_S03_L053_B1_T1_P022.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_94AF19F0_1F2A_41EC_8CB6_96C76227811F_X013_S03_L054_B1_T2_P027.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_0142AAAC_FFE8_43B7_AB99_02F7A1740567_X022_S06_L057_B1_T1_P050.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_0142AAAC_FFE8_43B7_AB99_02F7A1740567_X022_S06_L058_B1_T2_P051.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_257C1DDD_D646_4F59_BCEE_1403A58BBD80_X023_S06_L083_B1_T1_P056.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_257C1DDD_D646_4F59_BCEE_1403A58BBD80_X023_S06_L084_B1_T2_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_294C77E4_79D8_4383_8780_324F3DBB1839_X016_S02_L049_B1_T1_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_294C77E4_79D8_4383_8780_324F3DBB1839_X016_S02_L050_B1_T2_P035.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_359C8EB0_48BC_47D5_AD78_38B0AE6C5FB8_X018_S05_L032_B1_T1_P042.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_359C8EB0_48BC_47D5_AD78_38B0AE6C5FB8_X018_S05_L033_B1_T2_P043.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_474D1E91_69F9_4779_B8CE_18F4946A0D9A_X021_S03_L028_B1_T1_P054.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_474D1E91_69F9_4779_B8CE_18F4946A0D9A_X021_S03_L029_B1_T2_P054.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_852B450B_EBCC_4DD0_BB5A_FBF238C7BD89_X019_S02_L050_B1_T1_P041.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_852B450B_EBCC_4DD0_BB5A_FBF238C7BD89_X019_S02_L051_B1_T2_P044.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_873F57F9_5F1E_4B78_AE1A_65C9C6B487CE_X022_S07_L059_B1_T1_P050.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_873F57F9_5F1E_4B78_AE1A_65C9C6B487CE_X022_S07_L060_B1_T2_P047.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_1253F1AD_60FD_4536_97A8_E84B756E5E52_X004_S05_L051_B1_T1_P007.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_2106DD3A_1909_46E8_82BD_865333A12F8E_X013_S09_L065_B1_T1_P029.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_3403FB04_ED81_489B_991C_715B423B87DB_X014_S02_L003_B1_T1_P030.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_3403FB04_ED81_489B_991C_715B423B87DB_X014_S02_L004_B1_T2_P037.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_6108FBB3_DFCE_4A67_A16D_0547827F058C_X008_S04_L007_B1_T1_P018.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_6108FBB3_DFCE_4A67_A16D_0547827F058C_X008_S04_L008_B1_T2_P020.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_08499A64_3FD8_4E62_AF08_3C66AF93CAE7_X003_S05_L027_B1_T1_P003.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_08499A64_3FD8_4E62_AF08_3C66AF93CAE7_X003_S05_L028_B1_T2_P003.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8762FC8D_91FC_42B8_B74D_141298970EFC_X001_S03_L005_B1_T1_P004.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_8762FC8D_91FC_42B8_B74D_141298970EFC_X001_S03_L006_B1_T2_P004.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_08891EB8_2C41_400C_8A66_BEE040CCA71D_X011_S12_L023_B1_T1_P027.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_08891EB8_2C41_400C_8A66_BEE040CCA71D_X011_S12_L024_B1_T2_P029.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_037238B0_8FB6_4ECC_9970_93E84F9286EF_X005_S09_L020_B1_T1_P015.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_037238B0_8FB6_4ECC_9970_93E84F9286EF_X005_S09_L021_B1_T2_P014.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_73316C1E_DDE9_4C39_8AFD_5E8A1090C425_X006_S11_L052_B1_T1_P009.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_73316C1E_DDE9_4C39_8AFD_5E8A1090C425_X006_S11_L053_B1_T2_P010.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_503860D0_4069_4F6F_91F2_FFB24C3EA1E5_X005_S08_L018_B1_T1_P014.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_503860D0_4069_4F6F_91F2_FFB24C3EA1E5_X005_S08_L019_B1_T2_P015.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_01112370_4F6F_4A20_9BE0_7975C3465268_X017_S04_L007_B1_T1_P042.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_01112370_4F6F_4A20_9BE0_7975C3465268_X017_S04_L008_B1_T2_P044.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_3899957C_47D6_4274_AB3C_C9E6DE56437F_X008_S09_L017_B1_T1_P029.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_3899957C_47D6_4274_AB3C_C9E6DE56437F_X008_S09_L018_B1_T2_P021.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_94253529_096B_410E_88D2_80AA79361AD9_X015_S11_L044_B1_T1_P035.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A3A3CD8F_3497_4615_B288_14E6290DC441_X007_S02_L058_B1_T1_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A3A3CD8F_3497_4615_B288_14E6290DC441_X007_S02_L059_B1_T2_P010.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A8BD9D89_6D7F_4423_91A6_1627E36887F6_X022_S10_L065_B1_T1_P047.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A91AADEA_8299_46D9_A250_76896D690AFD_X016_S12_L069_B1_T1_P035.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A91AADEA_8299_46D9_A250_76896D690AFD_X016_S12_L070_B1_T2_P033.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A6462625_3970_44A4_B4CB_7B4AE86648CD_X013_S04_L055_B1_T1_P023.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_A6462625_3970_44A4_B4CB_7B4AE86648CD_X013_S04_L056_B1_T2_P025.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_BE0F871D_DFA0_4B91_A068_41A57C105E1B_X007_S05_L064_B1_T1_P013.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_BE061C50_495F_46AD_BF95_09C8AD8B5365_X022_S05_L056_B1_T1_P049.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_C9C8D426_A3FD_4455_89A9_768BC01D66A9_X009_S01_L025_B1_T1_P016.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_C9C8D426_A3FD_4455_89A9_768BC01D66A9_X009_S01_L026_B1_T2_P021.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_C147AAD5_A8F1_41D5_8709_21820BE50902_X016_S11_L067_B1_T1_P032.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_C147AAD5_A8F1_41D5_8709_21820BE50902_X016_S11_L068_B1_T2_P033.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CA5AB738_9366_4908_B573_92C041E15471_X020_S05_L009_B1_T1_P052.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CA5AB738_9366_4908_B573_92C041E15471_X020_S05_L010_B1_T2_P050.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CB96A542_7AC1_4FEC_A5D2_458D8EEDC6C4_X013_S05_L057_B1_T1_P023.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CB96A542_7AC1_4FEC_A5D2_458D8EEDC6C4_X013_S05_L058_B1_T2_P024.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CC102C17_C1CA_427A_8C7D_D3E79748A0CD_X004_S04_L049_B1_T1_P006.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CC102C17_C1CA_427A_8C7D_D3E79748A0CD_X004_S04_L050_B1_T2_P006.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CE1F0208_85C0_40DD_8C03_1595C8A6612C_X007_S01_L056_B1_T1_P031.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CE1F0208_85C0_40DD_8C03_1595C8A6612C_X007_S01_L057_B1_T2_P009.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CF183B78_9738_4D5F_8A57_82B8C8DD63AF_X009_S09_L041_B1_T1_P020.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_CF183B78_9738_4D5F_8A57_82B8C8DD63AF_X009_S09_L042_B1_T2_P018.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_D6F24369_2AB6_4A49_B3A4_B3B61D5A13B4_X014_S06_L011_B1_T1_P037.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_D6F24369_2AB6_4A49_B3A4_B3B61D5A13B4_X014_S06_L012_B1_T2_P066.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_D88D85E7_B4D4_4ACB_9815_996EE6B3650D_X011_S03_L005_B1_T1_P028.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_D88D85E7_B4D4_4ACB_9815_996EE6B3650D_X011_S03_L006_B1_T2_P029.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_DA8D04F8_C74F_4D50_B0B6_BBD15047FE72_X008_S02_L003_B1_T1_P018.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_DA8D04F8_C74F_4D50_B0B6_BBD15047FE72_X008_S02_L004_B1_T2_P020.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_DD69EDE9_142D_46E2_AA06_58D07D3230FB_X016_S08_L061_B1_T1_P035.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_DD69EDE9_142D_46E2_AA06_58D07D3230FB_X016_S08_L062_B1_T2_P034.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_DFFD80D2_34AB_4EE6_A7D4_06B557634CC4_X010_S07_L061_B1_T1_P017.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_DFFD80D2_34AB_4EE6_A7D4_06B557634CC4_X010_S07_L062_B1_T2_P016.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_E2534288_EFAE_4EF6_98A0_0A931BC20FF1_X017_S07_L013_B1_T1_P039.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_E2534288_EFAE_4EF6_98A0_0A931BC20FF1_X017_S07_L014_B1_T2_P040.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_EC8621C3_0A39_478F_BE88_2F1D1A19DA2D_X023_S07_L085_B1_T1_P056.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_EC8621C3_0A39_478F_BE88_2F1D1A19DA2D_X023_S07_L086_B1_T2_P055.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_EF17C882_9808_4676_9DFA_432D34290B33_X023_S15_L101_B1_T1_P056.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_EF17C882_9808_4676_9DFA_432D34290B33_X023_S15_L102_B1_T2_PMRG.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FB1C995E_6C78_414A_B74C_8C77CD924348_X015_S09_L041_B1_T1_P033.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FB1C995E_6C78_414A_B74C_8C77CD924348_X015_S09_L042_B1_T2_P032.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FB055B59_7512_40E4_8547_39798A4C9B8C_X011_S09_L017_B1_T1_P024.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FB055B59_7512_40E4_8547_39798A4C9B8C_X011_S09_L018_B1_T2_P025.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FCD2477B_7E05_4EB7_BD63_302496AEA537_X017_S11_L021_B1_T1_P041.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FCD2477B_7E05_4EB7_BD63_302496AEA537_X017_S11_L022_B1_T2_P040.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FE43880C_3F93_4463_9C91_5A2DE7130718_X009_S11_L045_B1_T1_P016.globalnorm_c0.0_l200_g30_peaks.narrowPeak",
        "/ifs/scratch/c2b2/ac_lab/jk3755/tcga/peaks/BRCA_FE43880C_3F93_4463_9C91_5A2DE7130718_X009_S11_L046_B1_T2_P018.globalnorm_c0.0_l200_g30_peaks.narrowPeak"
