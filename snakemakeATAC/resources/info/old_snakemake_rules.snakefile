
####
rule STEP4_bowtie2_refgenome_align:
     input:
         a="{path}preprocessing/2/{sample}_L{lane}_R1.good.fastq",
         b="{path}preprocessing/2/{sample}_L{lane}_R2.good.fastq"
     output:
         "{path}preprocessing/3/{sample}_L{lane}.sam"
     threads:
         15
     conda:
         "resources/envs/bowtie2.yaml"
     resources:
         mem_mb=lambda params, attempt: attempt * 25000,
         run_time=lambda params, attempt: attempt * 12
     shell:
         "bowtie2 -q -p {threads} -X 1000 -x genomes/bowtie2/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}_L{wildcards.lane}.align.txt"

####
rule STEP5_coordinate_sort_sam:
     input:
         "{path}preprocessing/3/{sample}_L{lane}.sam"
    output:
        "{path}preprocessing/4/{sample}_L{lane}.cs.sam"
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output} -O sam"

####
rule STEP12_mapq_filter:
    input:
        "{path}preprocessing/c/{sample}.dp.bam"
    output:
        "{path}preprocessing/d/{sample}.final.bam"
    conda:
        "resources/envs/samtools.yaml"
    resources:
        run_time=lambda params, attempt: attempt * 4
    shell:
        "samtools view -h -q 2 -b {input} > {output}"

####
rule MACS2_call_peaks_p01:
    input:
        a="{path}aligned/{sample}.bam",
        b="{path}aligned/{sample}.bam.bai"
    output:
        "{path}peaks/{sample}_p01_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}_p01 --outdir {wildcards.path}peaks/ --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
