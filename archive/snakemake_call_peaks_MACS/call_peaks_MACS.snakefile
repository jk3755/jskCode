configfile: "resources/config/config.yaml"

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

rule AGGREGATOR_peak_calling:
    input:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak"
    output:
        "{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.peaks"
    shell:
        "touch {output}"

rule MACS2_peaks_global_normilization_p001_one_sample:
    input:
        a="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}aligned/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/gn/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak"
    conda:
        "resources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm_p001 --outdir {wildcards.path}peaks/gn --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"
