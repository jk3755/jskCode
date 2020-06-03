
####
configfile: "config.yaml"

####
rule cosma_IDs:
    input:
        ancient(expand("data/{ID}/{ID}_1.fastq", ID=config["cosmaIDs"]))

####
rule fastq_dump_sra:
    output:
        "data/{SRA_ID}/{SRA_ID}_1.fastq"
    threads:
        20
    shell:
        "fasterq-dump -O data/{wildcards.SRA_ID} -e {threads} --split-files {wildcards.SRA_ID}"
