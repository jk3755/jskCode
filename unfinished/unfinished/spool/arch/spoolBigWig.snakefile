#############################################################################################################
#############################################################################################################
rule bigwig_cosma:
    input:
        "data/cosma/DAbaz2b/bigwig/DonorA_Baz2B.rep1.refhg38.bw",
        "data/cosma/DAluf/bigwig/DonorA_Luf.rep1.refhg38.bw",
        "data/cosma/DAprog/bigwig/DonorA_Progenitor.rep1.refhg38.bw",
        "data/cosma/DBbaz2b/bigwig/DonorB_Baz2B.rep1.refhg38.bw",
        "data/cosma/DBluf/bigwig/DonorB_Luf.rep1.refhg38.bw",
        "data/cosma/DBprog/bigwig/DonorB_Progenitor.rep1.refhg38.bw"
