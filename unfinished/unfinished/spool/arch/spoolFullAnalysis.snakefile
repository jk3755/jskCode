########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule full_lncap:
    input:
        "data/pros/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.full",
        "data/pros/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.full",
        "data/pros/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.full",
        "data/pros/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.full",
        "data/pros/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.full",
        "data/pros/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.full",
        "data/pros/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.full",
        "data/pros/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.full"


########################################################################################################################################
#### MISC #############################################################################################################################
########################################################################################################################################

rule full_temp:
    input:
        "data/pros/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.full",
        "data/pros/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.full",
        "data/pros/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.full",
        "data/pros/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.full",
        "data/pros/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.full",
        "data/pros/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.full",
        "data/pros/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.full",
        "data/pros/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.full",
        "data/cosma/ex01/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.preprocessing",
        "data/cosma/ex01/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.preprocessing",
        "data/cosma/ex01/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.preprocessing",
        "data/cosma/ex01/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.preprocessing",
        "data/cosma/ex01/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.preprocessing",
        "data/cosma/ex01/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.preprocessing"
