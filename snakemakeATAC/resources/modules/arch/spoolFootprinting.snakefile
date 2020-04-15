#############################################################################################################
#############################################################################################################
rule footprinting_lncap_pwm95:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_wt01_pwm95:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_wt02_pwm95:
    input:
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cr01_pwm95:
    input:
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cr02_pwm95:
    input:
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cr04_pwm95:
    input:
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cr05_pwm95:
    input:
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cr07_pwm95:
    input:
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cr08_pwm95:
    input:
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

#############################################################################################################
#############################################################################################################
rule wt01AR:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

rule cr01AR:
    input:
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

rule cr02AR:
    input:
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

rule cr04AR:
    input:
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

rule cr05AR:
    input:
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

rule cr07AR:
    input:
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

rule cr08AR:
    input:
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.pwm90.uncorrected.sm"

############################################################################################################
############################################################################################################
rule corrected_test:
    input:
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.pwm95.corrected.sm"





############################################################################################################
############################################################################################################
rule run_WGS:
    input:
        "data/wgs/bigwig/SRR1554090.bw",
        "data/wgs/operations/aggregators/SRR1554090.footprinting.pwm90.wgs.lncap.complete",
        "data/wgs/operations/aggregators/SRR1554090.footprinting.pwm90.wgs.cosma.complete"

rule bigwig_WGS:
    input:
        "data/wgs/bigwig/SRR1554090.bw"

rule lncap_WGS:
    input:
        "data/wgs/operations/aggregators/SRR1554090.footprinting.pwm90.wgs.lncap.complete"

rule cosma_WGS:
    input:
        "data/wgs/operations/aggregators/SRR1554090.footprinting.pwm90.wgs.cosma.complete"
