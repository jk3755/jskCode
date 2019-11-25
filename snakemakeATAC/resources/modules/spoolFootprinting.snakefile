########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule footprinting_lncap:
    input:
        "data/pros/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.uncorrected",
        "data/pros/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.uncorrected"

rule footprinting_lncap_chr:
    input:
        "data/pros/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.uncorrected.chr",
        "data/pros/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.uncorrected.chr"
