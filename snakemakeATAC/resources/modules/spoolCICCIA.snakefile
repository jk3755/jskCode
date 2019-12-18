#############################################################################################################
#############################################################################################################
rule ciccia_test:
    input:
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

"{path}operations/aggregators/{sample}.rep{repnum}.ref{refgenome}.footprinting.minimum{minsites}.chunks{totalchunk}.uncorrected.sm"
