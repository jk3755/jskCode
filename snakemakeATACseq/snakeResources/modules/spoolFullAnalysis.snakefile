## run preprocessing and footprinting analysis for all samples to date

rule full_ls1034:
    input:
        "data/coad/ls1034/wt04/operations/LNCaP-WT-04.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt05/operations/LNCaP-WT-05.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt06/operations/LNCaP-WT-06.rep1.refhg38.full_analysis.finished"

rule full_lncap:
    input:
        "data/pros/lncap/wt01/operations/LNCaP-WT-01.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/wt02/operations/LNCaP-WT-02.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr01/operations/LNCaP-CR-01.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr02/operations/LNCaP-CR-02.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr04/operations/LNCaP-CR-04.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr05/operations/LNCaP-CR-05.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr07/operations/LNCaP-CR-07.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr08/operations/LNCaP-CR-08.rep1.refhg38.full_analysis.finished"

rule full_lncap1:
    input:
        "data/pros/lncap/wt01/operations/LNCaP-WT-01.rep1.refhg38.full_analysis.finished"

rule full_lncap2:
    input:
        "data/pros/lncap/wt02/operations/LNCaP-WT-02.rep1.refhg38.full_analysis.finished"

rule full_lncap3:
    input:
        "data/pros/lncap/cr01/operations/LNCaP-CR-01.rep1.refhg38.full_analysis.finished"

rule full_lncap4:
    input:
        "data/pros/lncap/cr02/operations/LNCaP-CR-02.rep1.refhg38.full_analysis.finished"

rule full_lncap5:
    input:
        "data/pros/lncap/cr04/operations/LNCaP-CR-04.rep1.refhg38.full_analysis.finished"

rule full_lncap6:
    input:
        "data/pros/lncap/cr05/operations/LNCaP-CR-05.rep1.refhg38.full_analysis.finished"

rule full_lncap7:
    input:
        "data/pros/lncap/cr07/operations/LNCaP-CR-07.rep1.refhg38.full_analysis.finished"

rule full_lncap8:
    input:
        "data/pros/lncap/cr08/operations/LNCaP-CR-08.rep1.refhg38.full_analysis.finished"

rule lncap:
    input:
        "data/pros/lncap/wt01/operations/LNCaP-WT-01.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/wt02/operations/LNCaP-WT-02.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/cr01/operations/LNCaP-CR-01.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/cr02/operations/LNCaP-CR-02.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/cr04/operations/LNCaP-CR-04.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/cr05/operations/LNCaP-CR-05.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/cr07/operations/LNCaP-CR-07.rep1.refhg38.generate_figures_preprocessing.complete",
        "data/pros/lncap/cr08/operations/LNCaP-CR-08.rep1.refhg38.generate_figures_preprocessing.complete"
