##################################################################################################################################################
##################################################################################################################################################
## Make sure to change target genes in config file
rule ciccia_ctrl_untr_r1:
    input:
        "data/ciccia/CTRL_UNTR_R1/operations/aggregators/MDA_MB436_CTRL_UNTR.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule ciccia_ctrl_untr_r2:
    input:
        "data/ciccia/CTRL_UNTR_R2/operations/aggregators/MDA_MB436_CTRL_UNTR.rep2.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule ciccia_ctrl_untr_r3:
    input:
        "data/ciccia/CTRL_UNTR_R3/operations/aggregators/MDA_MB436_CTRL_UNTR.rep3.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule ciccia_ko_untr_r1:
    input:
        "data/ciccia/KO_UNTR_R1/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep1.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule ciccia_ko_untr_r2:
    input:
        "data/ciccia/KO_UNTR_R2/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep2.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"

rule ciccia_ko_untr_r3:
    input:
        "data/ciccia/KO_UNTR_R3/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep3.refhg38.footprinting.minimum50000.chunks10.uncorrected.sm"


##################################################################################################################################################
##################################################################################################################################################
## Make sure to change target genes in config file
rule ciccia_ctrl_untr_r1_pwm80:
    input:
        "data/ciccia/CTRL_UNTR_R1/operations/aggregators/MDA_MB436_CTRL_UNTR.rep1.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"

rule ciccia_ctrl_untr_r2_pwm80:
    input:
        "data/ciccia/CTRL_UNTR_R2/operations/aggregators/MDA_MB436_CTRL_UNTR.rep2.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"

rule ciccia_ctrl_untr_r3_pwm80:
    input:
        "data/ciccia/CTRL_UNTR_R3/operations/aggregators/MDA_MB436_CTRL_UNTR.rep3.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"

rule ciccia_ko_untr_r1_pwm80:
    input:
        "data/ciccia/KO_UNTR_R1/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep1.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"

rule ciccia_ko_untr_r2_pwm80:
    input:
        "data/ciccia/KO_UNTR_R2/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep2.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"

rule ciccia_ko_untr_r3_pwm80:
    input:
        "data/ciccia/KO_UNTR_R3/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep3.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"

##################################################################################################################################################
##################################################################################################################################################
## no peaks
## Make sure to change target genes in config file
rule ciccia_ctrl_untr_r1_pwm90_nopeak:
    input:
        "data/ciccia/CTRL_UNTR_R1/operations/aggregators/MDA_MB436_CTRL_UNTR.rep1.refhg38.footprinting.pwm90.nopeak"

rule ciccia_ctrl_untr_r2_pwm90_nopeak:
    input:
        "data/ciccia/CTRL_UNTR_R2/operations/aggregators/MDA_MB436_CTRL_UNTR.rep2.refhg38.footprinting.pwm90.nopeak"

rule ciccia_ctrl_untr_r3_pwm90_nopeak:
    input:
        "data/ciccia/CTRL_UNTR_R3/operations/aggregators/MDA_MB436_CTRL_UNTR.rep3.refhg38.footprinting.pwm90.nopeak"

rule ciccia_ko_untr_r1_pwm90_nopeak:
    input:
        "data/ciccia/KO_UNTR_R1/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep1.refhg38.footprinting.pwm90.nopeak"

rule ciccia_ko_untr_r2_pwm90_nopeak:
    input:
        "data/ciccia/KO_UNTR_R2/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep2.refhg38.footprinting.pwm90.nopeak"

rule ciccia_ko_untr_r3_pwm90_nopeak:
    input:
        "data/ciccia/KO_UNTR_R3/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep3.refhg38.footprinting.pwm90.nopeak"

####
rule ciccia_ap1:
    input:
        "data/ciccia/CTRL_UNTR_R1/operations/aggregators/MDA_MB436_CTRL_UNTR.rep1.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10",
        "data/ciccia/CTRL_UNTR_R2/operations/aggregators/MDA_MB436_CTRL_UNTR.rep2.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10",
        "data/ciccia/CTRL_UNTR_R3/operations/aggregators/MDA_MB436_CTRL_UNTR.rep3.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10",
        "data/ciccia/KO_UNTR_R1/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep1.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10",
        "data/ciccia/KO_UNTR_R2/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep2.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10",
        "data/ciccia/KO_UNTR_R3/operations/aggregators/MDA_MB436_SL1KO10_UNTR.rep3.refhg38.footprinting.pwm80.uncorrected.sm.totalchunk10"
