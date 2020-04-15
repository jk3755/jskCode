
rule STEP15_make_bigwig:
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
	