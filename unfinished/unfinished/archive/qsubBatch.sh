#!/bin/bash
#
qsub resources/config/qsub_wgs_lncap_fork1.sh
#sleep 60
qsub resources/config/qsubSubmit_WT01_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_WT02_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_CR01_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_CR02_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_CR04_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_CR05_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_CR07_95_allGenes_Grouped.sh
#sleep 60
qsub resources/config/qsubSubmit_CR08_95_allGenes_Grouped.sh
#sleep 60