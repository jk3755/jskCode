rule hahn2019:
    input:
        "data/SRR10153177/operations/SRR10153177.preprocessing_complete",
        "data/SRR10153179/operations/SRR10153179.preprocessing_complete"

rule test:
    input:
        "data/SRR10153177/operations/SRR10153177.preprocessing_complete"
