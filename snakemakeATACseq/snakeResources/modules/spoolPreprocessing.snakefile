rule preprocessing_ls1034:
    input:
        "data/coad/ls1034/wt04/operations/LS1034-WT-04.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt05/operations/LS1034-WT-05.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt06/operations/LS1034-WT-06.rep1.refhg38.preprocessing.complete"


rule preprocessing_coad_only:
    input:
        "data/coad/h508/wt04/operations/H508-WT-04.rep1.refhg38.preprocessing.complete",
        "data/coad/hs675t/wt01/operations/Hs675T-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt04/operations/LS1034-WT-04.rep1.refhg38.preprocessing.complete",
        "data/coad/mdst8/wt01/operations/MDST8-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/snu16/wt01/operations/SNU16-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/snu61/wt01/operations/SNU61-WT-01.rep1.refhg38.preprocessing.complete"

rule preprocessing_panc_only:
    input:
        "data/panc/capan1/split/wt01r1/operations/CAPANI-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/hpafii/split/wt01r1/operations/HPAFII-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/kp4/split/wt01r1/operations/KP4-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc1/split/wt01r1/operations/PANC1-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc0403/split/wt01r1/operations/PANC0403-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/patu8ss89/split/wt01r1/operations/PATU8SS89-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/pk45h/split/wt01r1/operations/PK45H-WT-01-RUN1.rep1.refhg38.preprocessing.complete"

rule preprocessing_coad:
    input:
        #"data/coad/h508/wt01/operations/H508-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/h508/wt02/operations/H508-WT-02.rep1.refhg38.preprocessing.complete",
        #"data/coad/h508/wt03/operations/H508-WT-03.rep1.refhg38.preprocessing.complete",
        "data/coad/h508/wt04/operations/H508-WT-04.rep1.refhg38.preprocessing.complete",
        #"data/coad/h508/wt05/operations/H508-WT-05.rep1.refhg38.preprocessing.complete",
        #"data/coad/h508/wt06/operations/H508-WT-06.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/hs675t/wt01/operations/Hs675T-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/hs675t/wt02/operations/Hs675T-WT-02.rep1.refhg38.preprocessing.complete",
        #"data/coad/hs675t/wt03/operations/Hs675T-WT-03.rep1.refhg38.preprocessing.complete",
        #
        #"data/coad/ls1034/wt01/operations/LS1034-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/ls1034/wt02/operations/LS1034-WT-02.rep1.refhg38.preprocessing.complete",
        #"data/coad/ls1034/wt03/operations/LS1034-WT-03.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt04/operations/LS1034-WT-04.rep1.refhg38.preprocessing.complete",
        #"data/coad/ls1034/wt05/operations/LS1034-WT-05.rep1.refhg38.preprocessing.complete",
        #"data/coad/ls1034/wt06/operations/LS1034-WT-06.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/mdst8/wt01/operations/MDST8-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/mdst8/wt02/operations/MDST8-WT-02.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/snu16/wt01/operations/SNU16-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/snu16/wt02/operations/SNU16-WT-01.rep2.refhg38.preprocessing.complete",
        #"data/coad/snu16/wt03/operations/SNU16-WT-01.rep3.refhg38.preprocessing.complete",
        #"data/coad/snu16/wt04/operations/SNU16-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/snu16/wt05/operations/SNU16-WT-01.rep2.refhg38.preprocessing.complete",
        #"data/coad/snu16/wt06/operations/SNU16-WT-01.rep3.refhg38.preprocessing.complete",
        #
        "data/coad/snu61/wt01/operations/SNU61-WT-01.rep1.refhg38.preprocessing.complete",
        #"data/coad/snu61/wt02/operations/SNU61-WT-01.rep2.refhg38.preprocessing.complete",
        #"data/coad/snu61/wt03/operations/SNU61-WT-01.rep3.refhg38.preprocessing.complete"
        
rule preprocessing_all:
    input:
        #### PROS ########################################################################
        "data/pros/lncap/wt01/operations/LNCaP-WT-01.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/wt02/operations/LNCaP-WT-02.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr01/operations/LNCaP-CR-01.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr02/operations/LNCaP-CR-02.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr04/operations/LNCaP-CR-04.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr05/operations/LNCaP-CR-05.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr07/operations/LNCaP-CR-07.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr08/operations/LNCaP-CR-08.rep1.refhg38.preprocessing.complete",
        #### COAD ########################################################################
        "data/coad/h508/wt01/operations/H508-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/h508/wt02/operations/H508-WT-02.rep1.refhg38.preprocessing.complete",
        "data/coad/h508/wt03/operations/H508-WT-03.rep1.refhg38.preprocessing.complete",
        "data/coad/h508/wt04/operations/H508-WT-04.rep1.refhg38.preprocessing.complete",
        "data/coad/h508/wt05/operations/H508-WT-05.rep1.refhg38.preprocessing.complete",
        "data/coad/h508/wt06/operations/H508-WT-06.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/hs675t/wt01/operations/Hs675T-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/hs675t/wt02/operations/Hs675T-WT-02.rep1.refhg38.preprocessing.complete",
        "data/coad/hs675t/wt03/operations/Hs675T-WT-03.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/ls1034/wt01/operations/LS1034-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt02/operations/LS1034-WT-02.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt03/operations/LS1034-WT-03.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt04/operations/LS1034-WT-04.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt05/operations/LS1034-WT-05.rep1.refhg38.preprocessing.complete",
        "data/coad/ls1034/wt06/operations/LS1034-WT-06.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/mdst8/wt01/operations/MDST8-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/mdst8/wt02/operations/MDST8-WT-02.rep1.refhg38.preprocessing.complete",
        #
        "data/coad/snu16/wt01/operations/SNU16-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/snu16/wt02/operations/SNU16-WT-01.rep2.refhg38.preprocessing.complete",
        "data/coad/snu16/wt03/operations/SNU16-WT-01.rep3.refhg38.preprocessing.complete",
        "data/coad/snu16/wt04/operations/SNU16-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/snu16/wt05/operations/SNU16-WT-01.rep2.refhg38.preprocessing.complete",
        "data/coad/snu16/wt06/operations/SNU16-WT-01.rep3.refhg38.preprocessing.complete",
        #
        "data/coad/snu61/wt01/operations/SNU61-WT-01.rep1.refhg38.preprocessing.complete",
        "data/coad/snu61/wt02/operations/SNU61-WT-01.rep2.refhg38.preprocessing.complete",
        "data/coad/snu61/wt03/operations/SNU61-WT-01.rep3.refhg38.preprocessing.complete",
        #### PANC ################################################################################
        "data/panc/capan1/split/wt01r1/operations/CAPANI-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/capan1/split/wt01r2/operations/CAPANI-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/capan1/split/wt02r1/operations/CAPANI-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/capan1/split/wt02r2/operations/CAPANI-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/capan1/split/wt03r1/operations/CAPANI-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/capan1/split/wt03r2/operations/CAPANI-WT-03-RUN2.rep1.refhg38.preprocessing.complete",
        #
        "data/panc/hpafii/split/wt01r1/operations/HPAFII-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/hpafii/split/wt01r2/operations/HPAFII-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/hpafii/split/wt02r1/operations/HPAFII-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/hpafii/split/wt02r2/operations/HPAFII-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/hpafii/split/wt03r1/operations/HPAFII-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/hpafii/split/wt03r2/operations/HPAFII-WT-03-RUN2.rep1.refhg38.preprocessing.complete",
        #
        "data/panc/kp4/split/wt01r1/operations/KP4-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/kp4/split/wt01r2/operations/KP4-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/kp4/split/wt02r1/operations/KP4-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/kp4/split/wt02r2/operations/KP4-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/kp4/split/wt03r1/operations/KP4-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/kp4/split/wt03r2/operations/KP4-WT-03-RUN2.rep1.refhg38.preprocessing.complete",
        #
        "data/panc/panc1/split/wt01r1/operations/PANC1-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc1/split/wt01r2/operations/PANC1-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/panc1/split/wt02r1/operations/PANC1-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc1/split/wt02r2/operations/PANC1-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/panc1/split/wt03r1/operations/PANC1-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc1/split/wt03r2/operations/PANC1-WT-03-RUN2.rep1.refhg38.preprocessing.complete",
        #
        "data/panc/panc0403/split/wt01r1/operations/PANC0403-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc0403/split/wt01r2/operations/PANC0403-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/panc0403/split/wt02r1/operations/PANC0403-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc0403/split/wt02r2/operations/PANC0403-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/panc0403/split/wt03r1/operations/PANC0403-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/panc0403/split/wt03r2/operations/PANC0403-WT-03-RUN2.rep1.refhg38.preprocessing.complete",
        #
        "data/panc/patu8ss89/split/wt01r1/operations/PATU8SS89-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/patu8ss89/split/wt01r2/operations/PATU8SS89-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/patu8ss89/split/wt02r1/operations/PATU8SS89-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/patu8ss89/split/wt02r2/operations/PATU8SS89-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/patu8ss89/split/wt03r1/operations/PATU8SS89-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/patu8ss89/split/wt03r2/operations/PATU8SS89-WT-03-RUN2.rep1.refhg38.preprocessing.complete",
        #
        "data/panc/pk45h/split/wt01r1/operations/PK45H-WT-01-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/pk45h/split/wt01r2/operations/PK45H-WT-01-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/pk45h/split/wt02r1/operations/PK45H-WT-02-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/pk45h/split/wt02r2/operations/PK45H-WT-02-RUN2.rep1.refhg38.preprocessing.complete",
        "data/panc/pk45h/split/wt03r1/operations/PK45H-WT-03-RUN1.rep1.refhg38.preprocessing.complete",
        "data/panc/pk45h/split/wt03r2/operations/PK45H-WT-03-RUN2.rep1.refhg38.preprocessing.complete"

rule preprocessing_lncap:
    input:
        "data/pros/lncap/wt01/operations/LNCaP-WT-01.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/wt02/operations/LNCaP-WT-02.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr01/operations/LNCaP-CR-01.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr02/operations/LNCaP-CR-02.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr04/operations/LNCaP-CR-04.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr05/operations/LNCaP-CR-05.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr07/operations/LNCaP-CR-07.rep1.refhg38.preprocessing.complete",
        "data/pros/lncap/cr08/operations/LNCaP-CR-08.rep1.refhg38.preprocessing.complete"