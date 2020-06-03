########################################################################################################################################
#### IMPORT MODULES ####################################################################################################################
########################################################################################################################################
configfile: "resources/config/config.yaml"
include: "resources/modules/spool_hahn2019.snakefile"

########################################################################################################################################
#### PREPROCESSING #####################################################################################################################
########################################################################################################################################

rule STEP1_build_directory_structure:
    output:
        "{path}operations/all_dirs.built"
    shell:
        """
        ####
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}aligned
        mkdir -p -v {wildcards.path}aligned_NF
        mkdir -p -v {wildcards.path}peaks_NF
        ####
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
        ####
        touch {output}
        """

rule STEP2_gunzip:
    input:
        a="{path}operations/all_dirs.built",
        b="{path}fastq/{sample}_R{read}.fastq.gz"        
    output:
        "{path}preprocessing/1/{sample}_R{read}.fastq"
    shell:
        "gunzip -c {input.b} > {output}"

rule STEP3_fastp_filtering:
    input:
        a="{path}preprocessing/1/{sample}_R1.fastq",
        b="{path}preprocessing/1/{sample}_R2.fastq"
    output:
        a="{path}preprocessing/2/{sample}_R1.good.fastq",
        b="{path}preprocessing/2/{sample}_R2.good.fastq"
    threads:
        5
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    conda:
    	"resources/envs/fastp.yaml"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.a} -O {output.b} -w {threads} -h {wildcards.path}metrics/{wildcards.sample}.quality.html -j {wildcards.path}metrics/{wildcards.sample}.quality.json"

rule STEP4_STAR_align:
    input:
        a="{path}preprocessing/2/{sample}_R1.good.fastq",
        b="{path}preprocessing/2/{sample}_R2.good.fastq"
    output:
        "{path}preprocessing/3/{sample}Aligned.out.sam"
    threads:
        20
    conda:
        "resources/envs/star.yaml"
    resources:
        mem_mb=50000
    shell:
        "STAR --genomeDir genomes/star/hg38 \
        --outFileNamePrefix {wildcards.path}preprocessing/3/{wildcards.sample} \
        --readFilesIn {input.a} {input.b} \
        --readFilesCommand cat \
        --runThreadN {threads} \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNoverReadLmax 0.1 \
        --alignIntronMax 1 \
        --alignMatesGapMax 1000 \
        --alignEndsType EndToEnd \
        --genomeLoad LoadAndRemove \
        --seedSearchStartLmax 30 \
        --outSAMattributes All \
        --outSAMunmapped None \
        --outSAMtype SAM \
        --outSAMheaderHD @HD VN:1.4"

rule STEP5_coordinate_sort_sam:
    input:
        "{path}preprocessing/3/{sample}Aligned.out.sam"
    output:
        "{path}preprocessing/4/{sample}.cs.sam"
    threads:
        10
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output} -O sam -@ {threads}"

rule STEP6_blacklist_filter_and_bam_conversion:
    input:
        "{path}preprocessing/4/{sample}.cs.sam"
    output:
        a="{path}preprocessing/5/{sample}.blacklist",
        b="{path}preprocessing/6/{sample}.blrm.bam"
    conda:
        "resources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/blacklist/hg38.blacklist.bed -U {output.b} -@ 4 {input}"

rule STEP7_chrM_contamination:
    input:
        "{path}preprocessing/6/{sample}.blrm.bam"
    output:
        a="{path}preprocessing/7/{sample}.mitochondrial",
        b="{path}preprocessing/8/{sample}.good.bam"
    conda:
        "resources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/mtdna/mtdna.extents.bed -U {output.b} -@ 4 {input}"

rule STEP8_add_rg_and_sort_bam:
    input:
        "{path}preprocessing/8/{sample}.good.bam"
    output:
        "{path}preprocessing/9/{sample}.tagged.bam"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID={wildcards.sample} \
        RGLB={wildcards.sample} \
        RGPL={wildcards.sample} \
        RGPU={wildcards.sample} \
        RGSM={wildcards.sample}"

rule STEP9_clean_bam:
    input:
        "{path}preprocessing/9/{sample}.tagged.bam"
    output:
        "{path}preprocessing/a/{sample}.clean.bam"
    conda:
        "resources/envs/picard.yaml"
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"

rule STEP10_remove_pcr_duplicates:
    input:
        "{path}preprocessing/a/{sample}.clean.bam"
    output:
        "{path}preprocessing/c/{sample}.dp.bam"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output} \
        M={wildcards.path}metrics/{wildcards.sample}.duplication.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

rule STEP11_move_bam:
    input:
        "{path}preprocessing/c/{sample}.dp.bam"
    output:
        "{path}aligned/{sample}.bam"
    shell:
        "cp {input} {output}"

rule STEP12_build_bam_index:
    input:
        ancient("{path}aligned/{sample}.bam")
    output:
        "{path}aligned/{sample}.bam.bai"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

rule STEP13_queryname_sort_sam:
    input:
        "{path}aligned/{sample}.bam"
    output:
        "{path}aligned/{sample}.nsort.bam"
    threads:
        10
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort -n {input} -o {output} -O bam -@ {threads}"

rule STEP14_GENRICH_call_peaks_q05:
    input:
        "{path}aligned/{sample}.nsort.bam"
    output:
        "{path}peaks/{sample}_q05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
    	"Genrich -t {input} -o {wildcards.path}peaks/{wildcards.sample}_q05.narrowPeak -j -q 0.05"

rule STEP14_GENRICH_call_peaks_p05:
    input:
        "{path}aligned/{sample}.nsort.bam"
    output:
        "{path}peaks/{sample}_p05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
    	"Genrich -t {input} -o {wildcards.path}peaks/{wildcards.sample}_p05.narrowPeak -j -p 0.05"

rule STEP15_deeptools_filter_NF:
    input:
        "{path}aligned/{sample}.bam"
    output:
        "{path}aligned_NF/{sample}_NF.bam",
        "{path}aligned_NF/{sample}_nucleosome_reads.bam"
    threads:
        20
    conda:
        "resources/envs/deeptools.yaml"
    shell:
        "alignmentSieve -b {input} -o {output} -p {threads} --filterMetrics {wildcards.path}/aligned_NF/NF_filter_metrics.txt --filteredOutReads {wildcards.path}/aligned_NF/{wildcards.sample}_nucleosome_reads.bam --verbose --minFragmentLength 0 --maxFragmentLength 147"

rule STEP16_build_bam_index_NF_reads:
    input:
        ancient("{path}aligned_NF/{sample}_NF.bam")
    output:
        "{path}aligned_NF/{sample}_NF.bam.bai"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

rule STEP17_build_bam_index_nucleosome_reads:
    input:
        ancient("{path}aligned_NF/{sample}_nucleosome_reads.bam")
    output:
        "{path}aligned_NF/{sample}_nucleosome_reads.bam.bai"
    conda:
        "resources/envs/picard.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

rule STEP18_queryname_sort_sam_NF:
    input:
        "{path}aligned_NF/{sample}_NF.bam"
    output:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    threads:
        10
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort -n {input} -o {output} -O bam -@ {threads}"

rule STEP19_queryname_sort_sam_nucleosome_reads:
    input:
        "{path}aligned_NF/{sample}_nucleosome_reads.bam"
    output:
        "{path}aligned_NF/{sample}_nucleosome_reads.nsort.bam"
    threads:
        10
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort -n {input} -o {output} -O bam -@ {threads}"

rule STEP20_GENRICH_call_peaks_q05_NF:
    input:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    output:
        "{path}peaks_NF/{sample}_NF_q05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
    	"Genrich -t {input} -o {wildcards.path}peaks_NF/{wildcards.sample}_NF_q05.narrowPeak -j -q 0.05"

rule STEP21_GENRICH_call_peaks_q05_nucleosome_free:
    input:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    output:
        "{path}peaks_NF/{sample}_nucleosome_reads_q05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
    	"Genrich -t {input} -o {wildcards.path}peaks_NF/{wildcards.sample}_nucleosome_reads_q05.narrowPeak -j -q 0.05"

rule STEP20_GENRICH_call_peaks_p05_NF:
    input:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    output:
        "{path}peaks_NF/{sample}_NF_p05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
    	"Genrich -t {input} -o {wildcards.path}peaks_NF/{wildcards.sample}_NF_p05.narrowPeak -j -p 0.05"

rule STEP21_GENRICH_call_peaks_p05_nucleosome_free:
    input:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    output:
        "{path}peaks_NF/{sample}_nucleosome_reads_p05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
    	"Genrich -t {input} -o {wildcards.path}peaks_NF/{wildcards.sample}_nucleosome_reads_p05.narrowPeak -j -p 0.05"

rule AGGREGATE_preprocessing:
    input:
        "{path}aligned/{sample}.bam.bai",
        "{path}aligned_NF/{sample}_NF.bam.bai",
        "{path}aligned_NF/{sample}_nucleosome_reads.bam.bai",
        "{path}peaks/{sample}_q05.narrowPeak",
        "{path}peaks_NF/{sample}_NF_q05.narrowPeak",
        "{path}peaks_NF/{sample}_nucleosome_reads_q05.narrowPeak",
        "{path}peaks/{sample}_p05.narrowPeak",
        "{path}peaks_NF/{sample}_NF_p05.narrowPeak"
    output:
        "{path}operations/{sample}.preprocessing_complete"
    shell:
        "touch {output}"