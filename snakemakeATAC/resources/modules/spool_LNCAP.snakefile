rule lncap1:
    input:
        "data/lncap/wt01/operations/LNCAP_WT01.preprocessing_complete",
        "data/lncap/wt02/operations/LNCAP_WT02.preprocessing_complete",
        "data/lncap/cr01/operations/LNCAP_CR01.preprocessing_complete",
        "data/lncap/cr02/operations/LNCAP_CR02.preprocessing_complete"

rule lncap2:
    input:
        "data/lncap/cr04/operations/LNCAP_CR04.preprocessing_complete",
        "data/lncap/cr05/operations/LNCAP_CR05.preprocessing_complete",
        "data/lncap/cr07/operations/LNCAP_CR07.preprocessing_complete",
        "data/lncap/cr08/operations/LNCAP_CR08.preprocessing_complete"