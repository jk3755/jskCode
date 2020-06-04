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


rule STEP15_deeptools_filter_NF:
    input:
        a="{path}aligned/{sample}.bam",
        b="{path}aligned/{sample}.bam.bai"
    output:
        a="{path}aligned_NF/{sample}_NF.bam",
        b="{path}aligned_NF/{sample}_nucleosome_reads.bam",
        c="{path}aligned_NF/{sample}_NF_filter_metrics.txt"
    threads:
        10
    conda:
        "resources/envs/deeptools.yaml"
    resources:
        mem_mb=50000
    shell:
        "alignmentSieve -b {input.a} -o {output.a} -p {threads} --filterMetrics {output.c} --filteredOutReads {output.b} --verbose --minFragmentLength 0 --maxFragmentLength 147"

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