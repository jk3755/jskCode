########################################################################################################################################
#### SPOOLING RULES ####################################################################################################################
########################################################################################################################################
rule human_SRR1554094:
    input:
        "data/human/bam/SRR1554094.bam.bai"

rule lncap_all:
    input:
        "data/lncap/bigwig/SRR5233717.bw",
        "data/lncap/bigwig/SRR5259501.bw",
        "data/lncap/bigwig/SRR5259502.bw"
        # 503 is not avail on sra
        #"data/lncap/bam/SRR5259503.bam.bai"

########################################################################################################################################
#### PROCESSING RULES ##################################################################################################################
########################################################################################################################################

## Build the directory structure
rule STEP1_build_dir:
    output:
        "{path}operations/main_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}goodfastq
        mkdir -p -v {wildcards.path}rawsam
        mkdir -p -v {wildcards.path}cssam
        mkdir -p -v {wildcards.path}rawbam
        mkdir -p -v {wildcards.path}bamrg
        mkdir -p -v {wildcards.path}cleanbam
        mkdir -p -v {wildcards.path}bammerged
        mkdir -p -v {wildcards.path}dedupbam
        mkdir -p -v {wildcards.path}bammapq
        mkdir -p -v {wildcards.path}bam
        mkdir -p -v {wildcards.path}bigwig
        touch {output}
        """

## Fastq QC filtering
rule STEP2_fastp_filtering:
    input:
    	"{path}operations/main_dir.built",
        a="{path}fastq/{sample}_1.part_{part}.fastq",
        b="{path}fastq/{sample}_2.part_{part}.fastq"
    output:
        c="{path}goodfastq/{sample}_1.part_{part}.gfq",
        d="{path}goodfastq/{sample}_2.part_{part}.gfq"
    threads:
        6
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    conda:
    	"resources/envs/fastp.yaml"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.c} -O {output.d} -w {threads}"

## Align reads to reference genome
rule STEP3_refgenome_align:
    input:
        a="{path}goodfastq/{sample}_1.part_{part}.gfq",
        b="{path}goodfastq/{sample}_2.part_{part}.gfq"
    output:
        "{path}rawsam/{sample}.part_{part}.sam"
    threads:
        6
    conda:
        "resources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 25000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "bowtie2 -q -p {threads} -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output}"

## Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP4_coordinate_sort_sam:
    input:
        "{path}rawsam/{sample}.part_{part}.sam"
    output:
        "{path}cssam/{sample}.part_{part}.sort.sam"
    threads:
        1
    conda:
        "resources/envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output} -O sam"
    
## Remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
rule STEP5_blacklist_filter_bam_conversion:
    input:
        "{path}cssam/{sample}.part_{part}.sort.sam"
    output:
        a="{path}rawbam/{sample}.part_{part}.blacklist.bam",
        b="{path}rawbam/{sample}.part_{part}.blrm.bam"
    conda:
        "resources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"

## Add @RG tags to the reads and perform coordinate sorting
rule STEP6_add_rg_and_cs_bam:
    input:
         "{path}rawbam/{sample}.part_{part}.blrm.bam"
    output:
        "{path}bamrg/{sample}.part_{part}.readgroup.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID={wildcards.sample} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU={wildcards.sample} \
        RGSM={wildcards.sample}"
    
## Clean the bam file
rule STEP7_clean_bam:
    input:
        "{path}bamrg/{sample}.part_{part}.readgroup.bam"
    output:
        "{path}cleanbam/{sample}.part_{part}.cleaned.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"

## Purge PCR duplicate reads. See below for merging rules, which are sample-specific
rule STEP9_purge_pcr_duplicates:
    input:
        "{path}bammerged/{sample}.merged.bam"
    output:
        "{path}dedupbam/{sample}.dp.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output} \
        M=duplication.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
    
## Filter reads for only uniquely mapping
rule STEP10_mapq_filter:
    input:
        "{path}dedupbam/{sample}.dp.bam"
    output:
        "{path}bammapq/{sample}.mqfiltered.bam"
    conda:
        "resources/envs/samtools.yaml"
    threads:
        1
    resources:
        run_time=lambda params, attempt: attempt * 4
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
    
## Move (copy) bam to parent directory and rename
rule STEP11_move_bam:
    input:
        "{path}bammapq/{sample}.mqfiltered.bam"
    output:
        "{path}bam/{sample}.bam"
    shell:
        "cp {input} {output}"

## Build the .bai index for the processed bam file
rule STEP12_build_bai_index:
    input:
        "{path}bam/{sample}.bam"
    output:
        "{path}bam/{sample}.bam.bai"
    conda:
        "resources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

## Make a bigwig file from the bam file
rule STEP13_make_bigwig:
    input:
        a="{path}bam/{sample}.bam",
        b="{path}bam/{sample}.bam.bai"
    output:
        "{path}bigwig/{sample}.bw"
    conda:
        "resources/envs/deeptools.yaml"
    threads:
        10
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p {threads} -v"

########################################################################################################################################
#### SAMPLE SPECIFIC MERGE RULES #######################################################################################################
########################################################################################################################################

#### human SRR1554094 ####
rule human_SRR1554094_merge_bam:
    input:
        p1="{path}cleanbam/SRR1554094.part_001.cleaned.bam",
        p2="{path}cleanbam/SRR1554094.part_002.cleaned.bam",
        p3="{path}cleanbam/SRR1554094.part_003.cleaned.bam",
        p4="{path}cleanbam/SRR1554094.part_004.cleaned.bam",
        p5="{path}cleanbam/SRR1554094.part_005.cleaned.bam",
        p6="{path}cleanbam/SRR1554094.part_006.cleaned.bam",
        p7="{path}cleanbam/SRR1554094.part_007.cleaned.bam",
        p8="{path}cleanbam/SRR1554094.part_008.cleaned.bam",
        p9="{path}cleanbam/SRR1554094.part_009.cleaned.bam",
        p10="{path}cleanbam/SRR1554094.part_010.cleaned.bam",
        p11="{path}cleanbam/SRR1554094.part_011.cleaned.bam",
        p12="{path}cleanbam/SRR1554094.part_012.cleaned.bam",
        p13="{path}cleanbam/SRR1554094.part_013.cleaned.bam",
        p14="{path}cleanbam/SRR1554094.part_014.cleaned.bam",
        p15="{path}cleanbam/SRR1554094.part_015.cleaned.bam",
        p16="{path}cleanbam/SRR1554094.part_016.cleaned.bam",
        p17="{path}cleanbam/SRR1554094.part_017.cleaned.bam",
        p18="{path}cleanbam/SRR1554094.part_018.cleaned.bam",
        p19="{path}cleanbam/SRR1554094.part_019.cleaned.bam",
        p20="{path}cleanbam/SRR1554094.part_020.cleaned.bam",
        p21="{path}cleanbam/SRR1554094.part_021.cleaned.bam",
        p22="{path}cleanbam/SRR1554094.part_022.cleaned.bam",
        p23="{path}cleanbam/SRR1554094.part_023.cleaned.bam",
        p24="{path}cleanbam/SRR1554094.part_024.cleaned.bam",
        p25="{path}cleanbam/SRR1554094.part_025.cleaned.bam",
        p26="{path}cleanbam/SRR1554094.part_026.cleaned.bam",
        p27="{path}cleanbam/SRR1554094.part_027.cleaned.bam",
        p28="{path}cleanbam/SRR1554094.part_028.cleaned.bam",
        p29="{path}cleanbam/SRR1554094.part_029.cleaned.bam",
        p30="{path}cleanbam/SRR1554094.part_030.cleaned.bam",
        p31="{path}cleanbam/SRR1554094.part_031.cleaned.bam",
        p32="{path}cleanbam/SRR1554094.part_032.cleaned.bam",
        p33="{path}cleanbam/SRR1554094.part_033.cleaned.bam",
        p34="{path}cleanbam/SRR1554094.part_034.cleaned.bam",
        p35="{path}cleanbam/SRR1554094.part_035.cleaned.bam",
        p36="{path}cleanbam/SRR1554094.part_036.cleaned.bam",
        p37="{path}cleanbam/SRR1554094.part_037.cleaned.bam",
        p38="{path}cleanbam/SRR1554094.part_038.cleaned.bam",
        p39="{path}cleanbam/SRR1554094.part_039.cleaned.bam",
        p40="{path}cleanbam/SRR1554094.part_040.cleaned.bam",
        p41="{path}cleanbam/SRR1554094.part_041.cleaned.bam",
        p42="{path}cleanbam/SRR1554094.part_042.cleaned.bam",
        p43="{path}cleanbam/SRR1554094.part_043.cleaned.bam",
        p44="{path}cleanbam/SRR1554094.part_044.cleaned.bam",
        p45="{path}cleanbam/SRR1554094.part_045.cleaned.bam",
        p46="{path}cleanbam/SRR1554094.part_046.cleaned.bam",
        p47="{path}cleanbam/SRR1554094.part_047.cleaned.bam",
        p48="{path}cleanbam/SRR1554094.part_048.cleaned.bam",
        p49="{path}cleanbam/SRR1554094.part_049.cleaned.bam",
        p50="{path}cleanbam/SRR1554094.part_050.cleaned.bam",
        p51="{path}cleanbam/SRR1554094.part_051.cleaned.bam",
        p52="{path}cleanbam/SRR1554094.part_052.cleaned.bam",
        p53="{path}cleanbam/SRR1554094.part_053.cleaned.bam",
        p54="{path}cleanbam/SRR1554094.part_054.cleaned.bam",
        p55="{path}cleanbam/SRR1554094.part_055.cleaned.bam",
        p56="{path}cleanbam/SRR1554094.part_056.cleaned.bam",
        p57="{path}cleanbam/SRR1554094.part_057.cleaned.bam",
        p58="{path}cleanbam/SRR1554094.part_058.cleaned.bam",
        p59="{path}cleanbam/SRR1554094.part_059.cleaned.bam",
        p60="{path}cleanbam/SRR1554094.part_060.cleaned.bam",
        p61="{path}cleanbam/SRR1554094.part_061.cleaned.bam",
        p62="{path}cleanbam/SRR1554094.part_062.cleaned.bam",
        p63="{path}cleanbam/SRR1554094.part_063.cleaned.bam",
        p64="{path}cleanbam/SRR1554094.part_064.cleaned.bam",
        p65="{path}cleanbam/SRR1554094.part_065.cleaned.bam",
        p66="{path}cleanbam/SRR1554094.part_066.cleaned.bam",
        p67="{path}cleanbam/SRR1554094.part_067.cleaned.bam",
        p68="{path}cleanbam/SRR1554094.part_068.cleaned.bam",
        p69="{path}cleanbam/SRR1554094.part_069.cleaned.bam",
        p70="{path}cleanbam/SRR1554094.part_070.cleaned.bam",
        p71="{path}cleanbam/SRR1554094.part_071.cleaned.bam",
        p72="{path}cleanbam/SRR1554094.part_072.cleaned.bam",
        p73="{path}cleanbam/SRR1554094.part_073.cleaned.bam",
        p74="{path}cleanbam/SRR1554094.part_074.cleaned.bam",
        p75="{path}cleanbam/SRR1554094.part_075.cleaned.bam",
        p76="{path}cleanbam/SRR1554094.part_076.cleaned.bam",
        p77="{path}cleanbam/SRR1554094.part_077.cleaned.bam",
        p78="{path}cleanbam/SRR1554094.part_078.cleaned.bam",
        p79="{path}cleanbam/SRR1554094.part_079.cleaned.bam",
        p80="{path}cleanbam/SRR1554094.part_080.cleaned.bam",
        p81="{path}cleanbam/SRR1554094.part_081.cleaned.bam",
        p82="{path}cleanbam/SRR1554094.part_082.cleaned.bam",
        p83="{path}cleanbam/SRR1554094.part_083.cleaned.bam",
        p84="{path}cleanbam/SRR1554094.part_084.cleaned.bam",
        p85="{path}cleanbam/SRR1554094.part_085.cleaned.bam",
        p86="{path}cleanbam/SRR1554094.part_086.cleaned.bam",
        p87="{path}cleanbam/SRR1554094.part_087.cleaned.bam",
        p88="{path}cleanbam/SRR1554094.part_088.cleaned.bam",
        p89="{path}cleanbam/SRR1554094.part_089.cleaned.bam",
        p90="{path}cleanbam/SRR1554094.part_090.cleaned.bam",
        p91="{path}cleanbam/SRR1554094.part_091.cleaned.bam",
        p92="{path}cleanbam/SRR1554094.part_092.cleaned.bam",
        p93="{path}cleanbam/SRR1554094.part_093.cleaned.bam",
        p94="{path}cleanbam/SRR1554094.part_094.cleaned.bam",
        p95="{path}cleanbam/SRR1554094.part_095.cleaned.bam",
        p96="{path}cleanbam/SRR1554094.part_096.cleaned.bam",
        p97="{path}cleanbam/SRR1554094.part_097.cleaned.bam",
        p98="{path}cleanbam/SRR1554094.part_098.cleaned.bam",
        p99="{path}cleanbam/SRR1554094.part_099.cleaned.bam",
        p100="{path}cleanbam/SRR1554094.part_100.cleaned.bam"
    output:
        "{path}bammerged/SRR1554094.merged.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MergeSamFiles \
        I={input.p1} \
        I={input.p2} \
        I={input.p3} \
        I={input.p4} \
        I={input.p5} \
        I={input.p6} \
        I={input.p7} \
        I={input.p8} \
        I={input.p9} \
        I={input.p10} \
        I={input.p11} \
        I={input.p12} \
        I={input.p13} \
        I={input.p14} \
        I={input.p15} \
        I={input.p16} \
        I={input.p17} \
        I={input.p18} \
        I={input.p19} \
        I={input.p20} \
        I={input.p21} \
        I={input.p22} \
        I={input.p23} \
        I={input.p24} \
        I={input.p25} \
        I={input.p26} \
        I={input.p27} \
        I={input.p28} \
        I={input.p29} \
        I={input.p30} \
        I={input.p31} \
        I={input.p32} \
        I={input.p33} \
        I={input.p34} \
        I={input.p35} \
        I={input.p36} \
        I={input.p37} \
        I={input.p38} \
        I={input.p39} \
        I={input.p40} \
        I={input.p41} \
        I={input.p42} \
        I={input.p43} \
        I={input.p44} \
        I={input.p45} \
        I={input.p46} \
        I={input.p47} \
        I={input.p48} \
        I={input.p49} \
        I={input.p50} \
        I={input.p51} \
        I={input.p52} \
        I={input.p53} \
        I={input.p54} \
        I={input.p55} \
        I={input.p56} \
        I={input.p57} \
        I={input.p58} \
        I={input.p59} \
        I={input.p60} \
        I={input.p61} \
        I={input.p62} \
        I={input.p63} \
        I={input.p64} \
        I={input.p65} \
        I={input.p66} \
        I={input.p67} \
        I={input.p68} \
        I={input.p69} \
        I={input.p70} \
        I={input.p71} \
        I={input.p72} \
        I={input.p73} \
        I={input.p74} \
        I={input.p75} \
        I={input.p76} \
        I={input.p77} \
        I={input.p78} \
        I={input.p79} \
        I={input.p80} \
        I={input.p81} \
        I={input.p82} \
        I={input.p83} \
        I={input.p84} \
        I={input.p85} \
        I={input.p86} \
        I={input.p87} \
        I={input.p88} \
        I={input.p89} \
        I={input.p90} \
        I={input.p91} \
        I={input.p92} \
        I={input.p93} \
        I={input.p94} \
        I={input.p95} \
        I={input.p96} \
        I={input.p97} \
        I={input.p98} \
        I={input.p99} \
        I={input.p100} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

#### lncap SRR5233717 ####
rule lncap_SRR5233717_merge_bam:
    input:
        p1="{path}cleanbam/SRR5233717.part_001.cleaned.bam",
        p2="{path}cleanbam/SRR5233717.part_002.cleaned.bam",
        p3="{path}cleanbam/SRR5233717.part_003.cleaned.bam",
        p4="{path}cleanbam/SRR5233717.part_004.cleaned.bam",
        p5="{path}cleanbam/SRR5233717.part_005.cleaned.bam",
        p6="{path}cleanbam/SRR5233717.part_006.cleaned.bam",
        p7="{path}cleanbam/SRR5233717.part_007.cleaned.bam",
        p8="{path}cleanbam/SRR5233717.part_008.cleaned.bam",
        p9="{path}cleanbam/SRR5233717.part_009.cleaned.bam",
        p10="{path}cleanbam/SRR5233717.part_010.cleaned.bam",
        p11="{path}cleanbam/SRR5233717.part_011.cleaned.bam",
        p12="{path}cleanbam/SRR5233717.part_012.cleaned.bam",
        p13="{path}cleanbam/SRR5233717.part_013.cleaned.bam",
        p14="{path}cleanbam/SRR5233717.part_014.cleaned.bam",
        p15="{path}cleanbam/SRR5233717.part_015.cleaned.bam",
        p16="{path}cleanbam/SRR5233717.part_016.cleaned.bam",
        p17="{path}cleanbam/SRR5233717.part_017.cleaned.bam",
        p18="{path}cleanbam/SRR5233717.part_018.cleaned.bam",
        p19="{path}cleanbam/SRR5233717.part_019.cleaned.bam",
        p20="{path}cleanbam/SRR5233717.part_020.cleaned.bam"
    output:
        "{path}bammerged/SRR5233717.merged.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MergeSamFiles \
        I={input.p1} \
        I={input.p2} \
        I={input.p3} \
        I={input.p4} \
        I={input.p5} \
        I={input.p6} \
        I={input.p7} \
        I={input.p8} \
        I={input.p9} \
        I={input.p10} \
        I={input.p11} \
        I={input.p12} \
        I={input.p13} \
        I={input.p14} \
        I={input.p15} \
        I={input.p16} \
        I={input.p17} \
        I={input.p18} \
        I={input.p19} \
        I={input.p20} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

#### lncap SRR5259501 ####
rule lncap_SRR5259501_merge_bam:
    input:
        p1="{path}cleanbam/SRR5259501.part_001.cleaned.bam",
        p2="{path}cleanbam/SRR5259501.part_002.cleaned.bam",
        p3="{path}cleanbam/SRR5259501.part_003.cleaned.bam",
        p4="{path}cleanbam/SRR5259501.part_004.cleaned.bam",
        p5="{path}cleanbam/SRR5259501.part_005.cleaned.bam",
        p6="{path}cleanbam/SRR5259501.part_006.cleaned.bam",
        p7="{path}cleanbam/SRR5259501.part_007.cleaned.bam",
        p8="{path}cleanbam/SRR5259501.part_008.cleaned.bam",
        p9="{path}cleanbam/SRR5259501.part_009.cleaned.bam",
        p10="{path}cleanbam/SRR5259501.part_010.cleaned.bam",
        p11="{path}cleanbam/SRR5259501.part_011.cleaned.bam",
        p12="{path}cleanbam/SRR5259501.part_012.cleaned.bam",
        p13="{path}cleanbam/SRR5259501.part_013.cleaned.bam",
        p14="{path}cleanbam/SRR5259501.part_014.cleaned.bam",
        p15="{path}cleanbam/SRR5259501.part_015.cleaned.bam",
        p16="{path}cleanbam/SRR5259501.part_016.cleaned.bam",
        p17="{path}cleanbam/SRR5259501.part_017.cleaned.bam",
        p18="{path}cleanbam/SRR5259501.part_018.cleaned.bam",
        p19="{path}cleanbam/SRR5259501.part_019.cleaned.bam",
        p20="{path}cleanbam/SRR5259501.part_020.cleaned.bam"
    output:
        "{path}bammerged/SRR5259501.merged.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MergeSamFiles \
        I={input.p1} \
        I={input.p2} \
        I={input.p3} \
        I={input.p4} \
        I={input.p5} \
        I={input.p6} \
        I={input.p7} \
        I={input.p8} \
        I={input.p9} \
        I={input.p10} \
        I={input.p11} \
        I={input.p12} \
        I={input.p13} \
        I={input.p14} \
        I={input.p15} \
        I={input.p16} \
        I={input.p17} \
        I={input.p18} \
        I={input.p19} \
        I={input.p20} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

#### lncap SRR5259502 ####
rule lncap_SRR5259502_merge_bam:
    input:
        p1="{path}cleanbam/SRR5259502.part_001.cleaned.bam",
        p2="{path}cleanbam/SRR5259502.part_002.cleaned.bam",
        p3="{path}cleanbam/SRR5259502.part_003.cleaned.bam",
        p4="{path}cleanbam/SRR5259502.part_004.cleaned.bam",
        p5="{path}cleanbam/SRR5259502.part_005.cleaned.bam",
        p6="{path}cleanbam/SRR5259502.part_006.cleaned.bam",
        p7="{path}cleanbam/SRR5259502.part_007.cleaned.bam",
        p8="{path}cleanbam/SRR5259502.part_008.cleaned.bam",
        p9="{path}cleanbam/SRR5259502.part_009.cleaned.bam",
        p10="{path}cleanbam/SRR5259502.part_010.cleaned.bam",
        p11="{path}cleanbam/SRR5259502.part_011.cleaned.bam",
        p12="{path}cleanbam/SRR5259502.part_012.cleaned.bam",
        p13="{path}cleanbam/SRR5259502.part_013.cleaned.bam",
        p14="{path}cleanbam/SRR5259502.part_014.cleaned.bam",
        p15="{path}cleanbam/SRR5259502.part_015.cleaned.bam",
        p16="{path}cleanbam/SRR5259502.part_016.cleaned.bam",
        p17="{path}cleanbam/SRR5259502.part_017.cleaned.bam",
        p18="{path}cleanbam/SRR5259502.part_018.cleaned.bam",
        p19="{path}cleanbam/SRR5259502.part_019.cleaned.bam",
        p20="{path}cleanbam/SRR5259502.part_020.cleaned.bam"
    output:
        "{path}bammerged/SRR5259502.merged.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MergeSamFiles \
        I={input.p1} \
        I={input.p2} \
        I={input.p3} \
        I={input.p4} \
        I={input.p5} \
        I={input.p6} \
        I={input.p7} \
        I={input.p8} \
        I={input.p9} \
        I={input.p10} \
        I={input.p11} \
        I={input.p12} \
        I={input.p13} \
        I={input.p14} \
        I={input.p15} \
        I={input.p16} \
        I={input.p17} \
        I={input.p18} \
        I={input.p19} \
        I={input.p20} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

#### lncap SRR5259503 ####
rule lncap_SRR5259503_merge_bam:
    input:
        p1="{path}cleanbam/SRR5259503.part_001.cleaned.bam",
        p2="{path}cleanbam/SRR5259503.part_002.cleaned.bam",
        p3="{path}cleanbam/SRR5259503.part_003.cleaned.bam",
        p4="{path}cleanbam/SRR5259503.part_004.cleaned.bam",
        p5="{path}cleanbam/SRR5259503.part_005.cleaned.bam",
        p6="{path}cleanbam/SRR5259503.part_006.cleaned.bam",
        p7="{path}cleanbam/SRR5259503.part_007.cleaned.bam",
        p8="{path}cleanbam/SRR5259503.part_008.cleaned.bam",
        p9="{path}cleanbam/SRR5259503.part_009.cleaned.bam",
        p10="{path}cleanbam/SRR5259503.part_010.cleaned.bam",
        p11="{path}cleanbam/SRR5259503.part_011.cleaned.bam",
        p12="{path}cleanbam/SRR5259503.part_012.cleaned.bam",
        p13="{path}cleanbam/SRR5259503.part_013.cleaned.bam",
        p14="{path}cleanbam/SRR5259503.part_014.cleaned.bam",
        p15="{path}cleanbam/SRR5259503.part_015.cleaned.bam",
        p16="{path}cleanbam/SRR5259503.part_016.cleaned.bam",
        p17="{path}cleanbam/SRR5259503.part_017.cleaned.bam",
        p18="{path}cleanbam/SRR5259503.part_018.cleaned.bam",
        p19="{path}cleanbam/SRR5259503.part_019.cleaned.bam",
        p20="{path}cleanbam/SRR5259503.part_020.cleaned.bam"
    output:
        "{path}bammerged/SRR5259503.merged.bam"
    conda:
        "resources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "picard MergeSamFiles \
        I={input.p1} \
        I={input.p2} \
        I={input.p3} \
        I={input.p4} \
        I={input.p5} \
        I={input.p6} \
        I={input.p7} \
        I={input.p8} \
        I={input.p9} \
        I={input.p10} \
        I={input.p11} \
        I={input.p12} \
        I={input.p13} \
        I={input.p14} \
        I={input.p15} \
        I={input.p16} \
        I={input.p17} \
        I={input.p18} \
        I={input.p19} \
        I={input.p20} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"
