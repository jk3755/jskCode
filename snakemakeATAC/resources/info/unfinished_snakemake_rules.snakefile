rule STEP16_NFreads_convert_bam_sam:
    input:
        "{path}aligned/{sample}.bam"
    output:
        "{path}preprocessing/e/{sample}_nf.sam"
    threads:
        1
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools view -h {input} > {output}"

rule STEP16_NFreads_filter_sam:
    input:
        "{path}preprocessing/e/{sample}_nf.sam"
    output:
        "{path}preprocessing/e/{sample}_nf_reads.sam"
    threads:
        1
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools view -h {input} |gawk '$9 < 100 && $9 > -100' >{output}"

rule STEP17_NFreads_convert_filtered_to_bam:
    input:
        "{path}preprocessing/e/{sample}_nf_reads.sam"
    output:
        "{path}preprocessing/f/{sample}_NF.bam"
    threads:
        1
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools view -b -h  {input} -o {output}"

rule STEP18_move_NF_bam:
    input:
        "{path}preprocessing/f/{sample}_NF.bam"
    output:
        "{path}aligned_NF/{sample}_NF.bam"
    shell:
        "cp {input} {output}"

rule STEP19_build_NF_bam_index:
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

rule STEP20_queryname_sort_NF_bam:
    input:
        "{path}aligned_NF/{sample}_NF.bam.bai"
    output:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    threads:
        10
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort -n {input} -o {output} -O bam -@ {threads}"

rule STEP21_GENRICH_call_NF_peaks_q05:
    input:
        "{path}aligned_NF/{sample}_NF.nsort.bam"
    output:
        "{path}peaks_NF/{sample}_NF_q05.narrowPeak"
    conda:
        "resources/envs/genrich.yaml"
    shell:
        "Genrich -t {input} -o {wildcards.path}peaks/{wildcards.sample}_NF_q05.narrowPeak -j -q 0.05"