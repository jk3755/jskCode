########################################################################################################################################
#### IMPORT SPOOLING MODULES ###########################################################################################################
########################################################################################################################################
include: "resources/modules/spool_COSMA_SRA.snakefile"

########################################################################################################################################
#### DIRECTORY STRUCTURE ###############################################################################################################
########################################################################################################################################
# Directory generation aggregator rule
rule AGGREGATOR_build_directory_structure:
    input:
        "{path}operations/directories/main_dir.built",
        "{path}operations/directories/operations_dir.built",
        "{path}operations/directories/preprocessing_dir.built",
        "{path}operations/directories/peaks_dir.built",
        "{path}operations/directories/metrics_dir.built" 
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
        mkdir -p -v {wildcards.path}aligned
        mkdir -p -v {wildcards.path}bigwig
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
        touch {output}
        """

# For QC metrics files
rule DIR_metrics:
    output:
        "{path}operations/directories/metrics_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}metrics/fastq
        mkdir -p -v {wildcards.path}metrics/align
        mkdir -p -v {wildcards.path}metrics/duplication
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

########################################################################################################################################
#### PREPROCESSING #####################################################################################################################
########################################################################################################################################

####
rule STEP1_fastp_filtering:
    input:
        a="{path}fastq/{sample}_1.fastq",
        b="{path}fastq/{sample}_2.fastq",
        c="{path}operations/directories/all_dirs.built"
    output:
        a="{path}preprocessing/1/{sample}_1.good.fq",
        b="{path}preprocessing/1/{sample}_2.good.fq"
    threads:
        20
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    conda:
    	"resources/envs/fastp.yaml"
    benchmark:
        "{path}benchmark/{sample}.benchmark.fastp.txt"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.a} -O {output.b} -w {threads} -h {wildcards.path}metrics/fastq.quality.html -j {wildcards.path}metrics/fastq.quality.json"
  
####
# rule STEP2_refgenome_align:
#     input:
#         a=ancient("{path}preprocessing/1/{sample}_1.good.fq"),
#         b=ancient("{path}preprocessing/1/{sample}_2.good.fq")
#     output:
#         "{path}preprocessing/2/{sample}.sam"
#     threads:
#         20
#     conda:
#         "resources/envs/bowtie2.yaml"
#     benchmark:
#         "{path}benchmark/{sample}.benchmark.align.txt"
#     resources:
#         mem_mb=lambda params, attempt: attempt * 50000,
#         run_time=lambda params, attempt: attempt * 24
#     shell:
#         "bowtie2 -q -p {threads} -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/align/{wildcards.sample}.refhg38.al"

####
rule STEP3_refgenome_align:
    input:
        a="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        b="{path}preprocessing/2/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4/{sample}.rep{repnum}.ref{refgenome}_L{lane}.sam"
    threads:
        15
    conda:
        "resources/envs/star.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "STAR --genomeDir genomes/star/hg38 --outFileNamePrefix {path}preprocessing/4/{sample}.rep{repnum}.ref{refgenome}_L{lane} --readFilesIn {input.a} {input.b} --readFilesCommand cat --runThreadN {threads} --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --alignIntronMax 1 --alignEndsType EndToEnd --genomeLoad NoSharedMemory --seedSearchStartLmax 30 --limitBAMsortRAM 0 --outSAMattributes NH HI NM MD AS XS --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4"

####
rule STEP3_coordinate_sort_sam:
    input:
        ancient("{path}preprocessing/2/{sample}.sam")
    output:
        "{path}preprocessing/3/{sample}.cs.sam"
    conda:
        "resources/envs/samtools.yaml"
    benchmark:
        "{path}benchmark/{sample}.benchmark.coordinate_sort_sam.txt"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        """
        rm -rf {wildcards.path}preprocessing/3/*
        samtools sort {input} -o {output} -O sam
        """

####
rule STEP4_blacklist_filter_and_bam_conversion:
    input:
        ancient("{path}preprocessing/3/{sample}.cs.sam")
    output:
        a="{path}preprocessing/4/{sample}.blacklist",
        b="{path}preprocessing/5/{sample}.blrm.bam"
    conda:
        "resources/envs/samtools.yaml"
    benchmark:
        "{path}benchmark/{sample}.benchmark.blacklist_filter_bam_conversion.txt"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"
    
####
rule STEP5_chrM_contamination:
    input:
        ancient("{path}preprocessing/5/{sample}.blrm.bam")
    output:
        a="{path}preprocessing/6/{sample}.mitochondrial",
        b="{path}preprocessing/7/{sample}.good.bam"
    conda:
        "resources/envs/samtools.yaml"
    benchmark:
        "{path}benchmark/{sample}.benchmark.blacklist_filter_bam_conversion.txt"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/mtdna/mtdna.extents.bed -U {output.b} -@ 4 {input}"

####
rule STEP6_add_rg_and_sort_bam:
    input:
        ancient("{path}preprocessing/7/{sample}.good.bam")
    output:
        "{path}preprocessing/8/{sample}.tagged.bam"
    conda:
        "resources/envs/picard.yaml"
    benchmark:
        "{path}benchmark/{sample}.add_rg_tags_sort.txt"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.sample} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.sample}.{wildcards.sample} \
        RGSM={wildcards.sample}"
    
####
rule STEP7_clean_bam:
    input:
        ancient("{path}preprocessing/8/{sample}.tagged.bam")
    output:
        "{path}preprocessing/9/{sample}.clean.bam"
    conda:
        "resources/envs/picard.yaml"
    benchmark:
        "{path}benchmark/{sample}.clean_bam.txt"
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"

####
rule STEP8_remove_pcr_duplicates:
    input:
        ancient("{path}preprocessing/9/{sample}.clean.bam")
    output:
        "{path}preprocessing/a/{sample}.dp.bam"
    conda:
        "resources/envs/picard.yaml"
    benchmark:
        "{path}benchmark/{sample}.remove_pcr_duplicates.txt"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output} \
        M={wildcards.path}metrics/duplication/{wildcards.sample}.duplication \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

####
rule STEP9_move_bam:
    input:
        ancient("{path}preprocessing/a/{sample}.dp.bam")
    output:
        "{path}aligned/{sample}.bam"
    benchmark:
        "{path}benchmark/{sample}.move_bam.txt"
    shell:
        """
        cp {wildcards.path}preprocessing/a/*.bam {wildcards.path}aligned/{wildcards.sample}.bam
        touch {output}
        """

####
rule STEP10_build_bai_index:
    input:
        ancient("{path}aligned/{sample}.bam")
    output:
        "{path}aligned/{sample}.bam.bai"
    conda:
        "resources/envs/picard.yaml"
    benchmark:
        "{path}benchmark/{sample}.build_bai_index.txt"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"
    
####
rule STEP11_make_bigwig:
    input:
        a=ancient("{path}aligned/{sample}.bam"),
        b=ancient("{path}aligned/{sample}.bam.bai")
    output:
        "{path}bigwig/{sample}.bw"
    conda:
        "resources/envs/deeptools.yaml"
    benchmark:
        "{path}benchmark/{sample}.bigwig.txt"
    threads:
        10
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p {threads} -v"

####
rule STEP12_MACS2_peaks_global_normilization_p001:
    input:
        a=ancient("{path}aligned/{sample}.bam"),
        b=ancient("{path}aligned/{sample}.bam.bai")
    output:
        "{path}peaks/{sample}_globalnorm_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    benchmark:
        "{path}benchmark/{sample}.macs2_peaks.txt"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "macs2 callpeak -t {input.a} -f BAM -g hs -n {wildcards.sample}_globalnorm_p001 --outdir {wildcards.path}peaks/ --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"

####
rule AGGREGATOR_preprocessing:
    input:
        ancient("{path}operations/directories/all_dirs.built"),
        ancient("{path}peaks/{sample}_globalnorm_p001_peaks.narrowPeak"),
        ancient("{path}aligned/{sample}.bam.bai")
    output:
        "{path}operations/aggregators/{sample}.preprocessing_complete"
    shell:
        """
        touch {output}
        """

####
rule CLEANUP_preprocessing:
    input:
        ancient("{path}operations/aggregators/{sample}.preprocessing_complete")
    output:
        "{path}operations/{sample}.preprocessing_cleaned"
    shell:
        """
        rm -rf {wildcards.path}preprocessing
        touch {output}
        """