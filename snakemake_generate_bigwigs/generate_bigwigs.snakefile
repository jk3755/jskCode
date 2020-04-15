####
configfile: "config.yaml"

####
rule make_bigwig:
    input:
        a=ancient("data/{sample}/{sample}.bam"),
        b=ancient("data/{sample}/{sample}.bam.bai")
    output:
        "data/{sample}/{sample}_normCPM_bin{binsize}.bw"
    conda:
        "deeptools.yaml"
    threads:
        20
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "bamCoverage --bam {input.a} --outFileName {output} --outFileFormat bigwig --effectiveGenomeSize 2913022398 --normalizeUsing CPM --exactScaling --binSize {wildcards.binsize} --numberOfProcessors {threads}"

rule bin_aggregator:
    input:
        #ancient("data/{sample}/{sample}_normCPM_bin1.bw"),
        ancient("data/{sample}/{sample}_normCPM_bin10.bw")
        #ancient("data/{sample}/{sample}_normCPM_bin100.bw")
        #ancient("data/{sample}/{sample}_normCPM_bin1000.bw")
    output:
        "data/{sample}_normCPM_allbins.complete"
    shell:
    	"touch {output}"

rule sample_expander:
    input:
        ancient(expand("data/{sample}_normCPM_allbins.complete", sample=config["sampleID"]))