rule lncap:
    input:
        "data/wt01/operations/WT01.preprocessing_complete",
        "data/wt02/operations/WT02.preprocessing_complete",
        "data/cr01/operations/CR01.preprocessing_complete",
        "data/cr02/operations/CR02.preprocessing_complete",
        "data/cr04/operations/CR04.preprocessing_complete",
        "data/cr05/operations/CR05.preprocessing_complete",
        "data/cr07/operations/CR07.preprocessing_complete",
        "data/cr08/operations/CR08.preprocessing_complete"

rule test:
    input:
        "data/cr01/operations/CR01.preprocessing_complete" 
