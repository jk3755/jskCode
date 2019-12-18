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
rule footprinting_cosma_pwm95:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/cosma/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/cosma/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/cosma/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/cosma/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.footprinting.pwm95.uncorrected.sm",
        "data/cosma/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_cosma_test:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_DAbaz2b_pwm95:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_DAluf_pwm95:
    input:
        "data/cosma/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_DAprog_pwm95:
    input:
        "data/cosma/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_DBbaz2b_pwm95:
    input:
        "data/cosma/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_DBluf_pwm95:
    input:
        "data/cosma/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.footprinting.pwm95.uncorrected.sm"

rule footprinting_DBprog_pwm95:
    input:
        "data/cosma/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.footprinting.pwm95.uncorrected.sm"


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
