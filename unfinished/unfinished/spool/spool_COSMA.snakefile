##################################################################################################################################################
##################################################################################################################################################
## Make sure to change target genes in config file
rule cosma_donorA_baz2b:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.preprocessing"
        
rule cosma_donorA_luf:
    input:
        "data/cosma/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule cosma_donorA_prog:
    input:
        "data/cosma/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule cosma_donorB_baz2b:
    input:
        "data/cosma/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule cosma_donorB_luf:
    input:
        "data/cosma/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule cosma_donorB_prog:
    input:
        "data/cosma/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"
