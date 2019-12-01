########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule footprinting_lncap:
    input:
        "data/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected",
        "data/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.uncorrected",
        "data/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.uncorrected",
        "data/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.uncorrected",
        "data/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.uncorrected",
        "data/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.uncorrected",
        "data/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.uncorrected",
        "data/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.uncorrected"

rule footprinting_temp:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.uncorrected",
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.uncorrected",
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.jp",
        "data/cosma/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.footprinting.jp",
        "data/cosma/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.footprinting.jp",
        "data/cosma/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.footprinting.jp",
        "data/cosma/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.footprinting.jp",
        "data/cosma/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.footprinting.jp",
        "data/cosma/UNd/operations/aggregators/Undetermined.rep1.refhg38.footprinting.jp"

rule footprinting_lncap_forked:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.uncorrected.forked",
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.uncorrected.forked"

rule footprinting_cosma_forked:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.jp.forked",
        "data/cosma/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.footprinting.jp.forked",
        "data/cosma/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.footprinting.jp.forked",
        "data/cosma/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.footprinting.jp.forked",
        "data/cosma/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.footprinting.jp.forked",
        "data/cosma/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.footprinting.jp.forked",
        "data/cosma/UNd/operations/aggregators/Undetermined.rep1.refhg38.footprinting.jp.forked"

rule footprinting_lncap_chunked:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.uncorrected.chunked",
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.uncorrected.chunked"

rule footprinting_test:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.uncorrected.chunked"

rule footprinting_cosma_chunked_pwm95:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked",
        "data/cosma/DAluf/operations/aggregators/DonorA_Luf.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked",
        "data/cosma/DAprog/operations/aggregators/DonorA_Progenitor.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked",
        "data/cosma/DBbaz2b/operations/aggregators/DonorB_Baz2B.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked",
        "data/cosma/DBluf/operations/aggregators/DonorB_Luf.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked",
        "data/cosma/DBprog/operations/aggregators/DonorB_Progenitor.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked",
        "data/cosma/UNd/operations/aggregators/Undetermined.rep1.refhg38.footprinting.pwm95.jp.uncorrected.chunked"

rule footprinting_cosma_chunked_test:
    input:
        "data/cosma/DAbaz2b/operations/aggregators/DonorA_Baz2B.rep1.refhg38.footprinting.jp.uncorrected.chunked"

rule footprinting_lncap_chunked_pwm95:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/wt02/operations/aggregators/LNCaP-WT-02.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/cr01/operations/aggregators/LNCaP-CR-01.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/cr02/operations/aggregators/LNCaP-CR-02.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/cr04/operations/aggregators/LNCaP-CR-04.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/cr05/operations/aggregators/LNCaP-CR-05.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/cr07/operations/aggregators/LNCaP-CR-07.rep1.refhg38.footprinting.pwm95.uncorrected.chunked",
        "data/lncap/cr08/operations/aggregators/LNCaP-CR-08.rep1.refhg38.footprinting.pwm95.uncorrected.chunked"

rule footprinting_wt01_pwm90:
    input:
        "data/lncap/wt01/operations/aggregators/LNCaP-WT-01.rep1.refhg38.footprinting.pwm90.uncorrected.chunked"
