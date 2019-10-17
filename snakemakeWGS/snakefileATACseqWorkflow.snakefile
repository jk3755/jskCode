########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################
configfile: "snakeResources/config/config.yaml"
include: "snakeResources/modules/spoolFullAnalysis.snakefile"

########################################################################################################################################
#### DIRECTORY STRUCTURE ###############################################################################################################
########################################################################################################################################

rule DIR_main:
    output:
        "{path}operations/main_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}benchmark
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}correlation
        mkdir -p -v {wildcards.path}bam
        mkdir -p -v {wildcards.path}bigwig
        mkdir -p -v {wildcards.path}figures
        mkdir -p -v {wildcards.path}seqbias
        touch {output}
        """

rule DIR_benchmark:
    output:
        "{path}operations/benchmark_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}benchmark/preprocessing
        mkdir -p -v {wildcards.path}benchmark/preprocessing/gunzip {wildcards.path}benchmark/preprocessing/fastp {wildcards.path}benchmark/preprocessing/mycoalign
        mkdir -p -v {wildcards.path}benchmark/preprocessing/align {wildcards.path}benchmark/preprocessing/coordsortsam {wildcards.path}benchmark/preprocessing/bamconversion
        mkdir -p -v {wildcards.path}benchmark/preprocessing/removemitochondrial {wildcards.path}benchmark/preprocessing/addRG {wildcards.path}benchmark/preprocessing/cleansam
        mkdir -p -v {wildcards.path}benchmark/preprocessing/mergelanes {wildcards.path}benchmark/preprocessing/purgeduplicates {wildcards.path}benchmark/preprocessing/mapqfilter
        mkdir -p -v {wildcards.path}benchmark/preprocessing/buildindex {wildcards.path}benchmark/preprocessing/bigwig {wildcards.path}benchmark/preprocessing/peaks
        mkdir -p -v {wildcards.path}benchmark/metrics
        mkdir -p -v {wildcards.path}benchmark/correlation
        mkdir -p -v {wildcards.path}benchmark/saturation
        mkdir -p -v {wildcards.path}benchmark/footprints
        mkdir -p -v {wildcards.path}benchmark/footprints/raw
        touch {output}
        """

rule DIR_metrics:
    output:
        "{path}operations/metrics_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}metrics/saturation
        mkdir -p -v {wildcards.path}metrics/fastq
        mkdir -p -v {wildcards.path}metrics/myco
        mkdir -p -v {wildcards.path}metrics/align
        mkdir -p -v {wildcards.path}metrics/genomecov
        mkdir -p -v {wildcards.path}metrics/totalreads
        mkdir -p -v {wildcards.path}metrics/fragsize
        mkdir -p -v {wildcards.path}metrics/peakannotation
        mkdir -p -v {wildcards.path}metrics/duplication
        touch {output}
        """

rule DIR_preprocessing:
    output:
        "{path}operations/preprocessing_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing/2fastq
        mkdir -p -v {wildcards.path}preprocessing/3goodfastq
        mkdir -p -v {wildcards.path}preprocessing/4mycoalign
        mkdir -p -v {wildcards.path}preprocessing/5align
        mkdir -p -v {wildcards.path}preprocessing/6raw
        mkdir -p -v {wildcards.path}preprocessing/6raw/mitochondrial {wildcards.path}preprocessing/6raw/blacklist {wildcards.path}preprocessing/6raw/nonblacklist
        mkdir -p -v {wildcards.path}preprocessing/7rgsort
        mkdir -p -v {wildcards.path}preprocessing/8merged
        mkdir -p -v {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique
        mkdir -p -v {wildcards.path}preprocessing/11bigwig
        mkdir -p -v {wildcards.path}preprocessing/12saturation
        touch {output}
        """

rule DIR_saturation:
    output:
        "{path}operations/saturation_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing/12saturation/downsampled {wildcards.path}preprocessing/12saturation/downsampled/raw
        mkdir -p -v {wildcards.path}preprocessing/12saturation/downsampled/sorted {wildcards.path}preprocessing/12saturation/downsampled/deduplicated
        mkdir -p -v {wildcards.path}preprocessing/12saturation/duplication
        mkdir -p -v {wildcards.path}preprocessing/12saturation/peaks
        mkdir -p -v {wildcards.path}preprocessing/12saturation/footprints
        touch {output}
        """

rule DIR_footprints:
    output:
        "{path}operations/footprints_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}footprints/sample_specific
        mkdir -p -v {wildcards.path}footprints/sample_specific/sites
        mkdir -p -v {wildcards.path}footprints/sample_specific/insertion_matrix
        mkdir -p -v {wildcards.path}footprints/sample_specific/raw
        mkdir -p -v {wildcards.path}footprints/sample_specific/aggregated
        mkdir -p -v {wildcards.path}footprints/sample_merged
        mkdir -p -v {wildcards.path}footprints/sample_merged/sites
        mkdir -p -v {wildcards.path}footprints/sample_merged/insertion_matrix
        mkdir -p -v {wildcards.path}footprints/sample_merged/raw
        mkdir -p -v {wildcards.path}footprints/sample_merged/aggregated
        touch {output}
        """

rule DIR_peaks:
    output:
        "{path}operations/peaks_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}peaks/localnorm {wildcards.path}peaks/globalnorm {wildcards.path}peaks/sample_merged
        touch {output}
        """

rule DIR_figures:
    output:
        "{path}operations/figures_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}figures/fragmentsizes
        mkdir -p -v {wildcards.path}figures/peakideogram
        mkdir -p -v {wildcards.path}figures/insertionprobability
        mkdir -p -v {wildcards.path}figures/motifalignedheatmap
        mkdir -p -v {wildcards.path}figures/seqbias
        touch {output}
        """

########################################################################################################################################
#### AGGREGATOR RULES ##################################################################################################################
########################################################################################################################################

## This rule determines what is run in the full analysis spooling option
rule AGGREGATOR_full_analysis:
    input:
        "{path}operations/{sample}.preprocessing.complete"
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.footprinting_analysis.complete"
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.seqbias_models.complete"
    output:
        "{path}operations/{sample}.full_analysis.finished"
    shell:
        "touch {output}"

## This rule determines what is run in the directory building step
rule AGGREGATOR_build_directory_structure:
    input:
        "{path}operations/main_dir.built",
        "{path}operations/benchmark_dir.built",
        "{path}operations/metrics_dir.built",
        "{path}operations/preprocessing_dir.built",
        "{path}operations/saturation_dir.built",
        "{path}operations/footprints_dir.built",
        "{path}operations/peaks_dir.built",
        "{path}operations/figures_dir.built"
    output:
        "{path}operations/all_dirs.built"
    shell:
        "touch {output}"

## This rule determines what is run in the preprocessing spooling option
rule AGGREGATOR_preprocessing:
    input:
        "{path}operations/all_dirs.built",
        "{path}bam/{sample}.bam.bai"
        #"{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw",
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.peak_calling.complete",
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing_metrics.complete",
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.generate_figures_preprocessing.complete"
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.saturation_analysis.complete"
    output:
        "{path}operations/{sample}.preprocessing_preclean.complete"
    priority:
        1
    shell:
        "touch {output}"

## Clean up all the intermediate files just before touching the final flag file
rule AGGREGATOR_preprocessing_clean_intermediate_data:
    input:
        "{path}operations/{sample}.preprocessing_preclean.complete"  
    output:
        "{path}operations/{sample}.preprocessing.complete"
    shell:
        """
        rm -rf {wildcards.path}preprocessing/2fastq
        rm -rf {wildcards.path}preprocessing/3goodfastq
        rm -rf {wildcards.path}preprocessing/4mycoalign
        rm -rf {wildcards.path}preprocessing/5align
        rm -rf {wildcards.path}preprocessing/6raw
        rm -rf {wildcards.path}preprocessing/7rgsort
        rm -rf {wildcards.path}preprocessing/8merged
        rm -rf {wildcards.path}preprocessing/9dedup
        rm -rf {wildcards.path}preprocessing/10unique
        rm -rf {wildcards.path}preprocessing/11bigwig
        rm -rf {wildcards.path}preprocessing/saturation
        touch {output}
        """

## This rule determines what is run for the footprinting analysis
rule AGGREGATOR_footprinting_analysis:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_specific_footprint_analysis.complete",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_merged_footprint_analysis.complete"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.footprinting_analysis.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################

## Gunzip the fastq files
rule STEP1_gunzip:
    input:
        a="{path}fastq/{sample}_{read}.part_{part}.fastq",
        b="{path}operations/all_dirs.built"
    output:
        c="{path}preprocessing/2fastq/{sample}_{read}.part_{part}.fastq"
    threads:
        1
    shell:
        "cp {input.a} {output.c}"

## Fastq QC filtering
rule STEP2_fastp_filtering:
    input:
        a="{path}preprocessing/2fastq/{sample}_1.part_{part}.fastq",
        b="{path}preprocessing/2fastq/{sample}_2.part_{part}.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}_1.part_{part}.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}_2.part_{part}.good.fq"
    threads:
        6
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    conda:
    	"snakeResources/envs/fastp.yaml"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.c} -O {output.d} -w 6"

## Align reads to reference genome
rule STEP4_refgenome_align:
    input:
        a="{path}preprocessing/3goodfastq/{sample}_1.part_{part}.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}_2.part_{part}.good.fq"
    output:
        "{path}preprocessing/5align/{sample}.part_{part}.sam"
    threads:
        6
    conda:
        "snakeResources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "bowtie2 -q -p 6 -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output}"

## Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP5_coordsort_sam:
    input:
        "{path}preprocessing/5align/{sample}.part_{part}.sam"
    output:
        "{path}preprocessing/5align/{sample}.part_{part}.cs.sam"
    threads:
        1
    conda:
        "snakeResources/envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output} -O sam"
    
## Remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
rule STEP6_blacklistfilter_bamconversion:
    input:
        "{path}preprocessing/5align/{sample}.part_{part}.cs.sam"
    output:
        a="{path}preprocessing/6raw/blacklist/{sample}.part_{part}.blacklist.bam",
        b="{path}preprocessing/6raw/nonblacklist/{sample}.part_{part}.blrm.bam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"
    
## Remove reads mapping to mitochondrial DNA
rule STEP7_chrM_contamination:
    input:
        "{path}preprocessing/6raw/nonblacklist/{sample}.part_{part}.blrm.bam"
    output:
        a="{path}preprocessing/6raw/mitochondrial/{sample}.part_{part}.mitochondrial.bam",
        b="{path}preprocessing/6raw/{sample}.part_{part}.goodbam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        5
    shell:
        "samtools view -b -h -o {output.a} -L genomes/mtdna/mtdna.extents.bed -U {output.b} -@ 4 {input}"

## Add @RG tags to the reads and perform coordinate sorting
rule STEP8_addrgandcsbam:
    input:
         b="{path}preprocessing/6raw/{sample}.part_{part}.goodbam"
    output:
        "{path}preprocessing/7rgsort/{sample}.part_{part}.rg.cs.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 1
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.sample} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.sample} \
        RGSM={wildcards.sample}"
    
## Clean the bam file
rule STEP9_cleansam:
    input:
        "{path}preprocessing/7rgsort/{sample}.part_{part}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}.part_{part}.clean.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"
    
## Merge reads from different NextSeq lanes
rule STEP10_mergelanes:
    input:
        p1="{path}preprocessing/7rgsort/{sample}.part_001.clean.bam",
        p2="{path}preprocessing/7rgsort/{sample}.part_002.clean.bam",
        p3="{path}preprocessing/7rgsort/{sample}.part_003.clean.bam",
        p4="{path}preprocessing/7rgsort/{sample}.part_004.clean.bam",
        p5="{path}preprocessing/7rgsort/{sample}.part_005.clean.bam",
        p6="{path}preprocessing/7rgsort/{sample}.part_006.clean.bam",
        p7="{path}preprocessing/7rgsort/{sample}.part_007.clean.bam",
        p8="{path}preprocessing/7rgsort/{sample}.part_008.clean.bam",
        p9="{path}preprocessing/7rgsort/{sample}.part_009.clean.bam",
        p10="{path}preprocessing/7rgsort/{sample}.part_010.clean.bam",
        p11="{path}preprocessing/7rgsort/{sample}.part_011.clean.bam",
        p12="{path}preprocessing/7rgsort/{sample}.part_012.clean.bam",
        p13="{path}preprocessing/7rgsort/{sample}.part_013.clean.bam",
        p14="{path}preprocessing/7rgsort/{sample}.part_014.clean.bam",
        p15="{path}preprocessing/7rgsort/{sample}.part_015.clean.bam",
        p16="{path}preprocessing/7rgsort/{sample}.part_016.clean.bam",
        p17="{path}preprocessing/7rgsort/{sample}.part_017.clean.bam",
        p18="{path}preprocessing/7rgsort/{sample}.part_018.clean.bam",
        p19="{path}preprocessing/7rgsort/{sample}.part_019.clean.bam",
        p20="{path}preprocessing/7rgsort/{sample}.part_020.clean.bam",
        p21="{path}preprocessing/7rgsort/{sample}.part_021.clean.bam",
        p22="{path}preprocessing/7rgsort/{sample}.part_022.clean.bam",
        p23="{path}preprocessing/7rgsort/{sample}.part_023.clean.bam",
        p24="{path}preprocessing/7rgsort/{sample}.part_024.clean.bam",
        p25="{path}preprocessing/7rgsort/{sample}.part_025.clean.bam",
        p26="{path}preprocessing/7rgsort/{sample}.part_026.clean.bam",
        p27="{path}preprocessing/7rgsort/{sample}.part_027.clean.bam",
        p28="{path}preprocessing/7rgsort/{sample}.part_028.clean.bam",
        p29="{path}preprocessing/7rgsort/{sample}.part_029.clean.bam",
        p30="{path}preprocessing/7rgsort/{sample}.part_030.clean.bam",
        p31="{path}preprocessing/7rgsort/{sample}.part_031.clean.bam",
        p32="{path}preprocessing/7rgsort/{sample}.part_032.clean.bam",
        p33="{path}preprocessing/7rgsort/{sample}.part_033.clean.bam",
        p34="{path}preprocessing/7rgsort/{sample}.part_034.clean.bam",
        p35="{path}preprocessing/7rgsort/{sample}.part_035.clean.bam",
        p36="{path}preprocessing/7rgsort/{sample}.part_036.clean.bam",
        p37="{path}preprocessing/7rgsort/{sample}.part_037.clean.bam",
        p38="{path}preprocessing/7rgsort/{sample}.part_038.clean.bam",
        p39="{path}preprocessing/7rgsort/{sample}.part_039.clean.bam",
        p40="{path}preprocessing/7rgsort/{sample}.part_040.clean.bam",
        p41="{path}preprocessing/7rgsort/{sample}.part_041.clean.bam",
        p42="{path}preprocessing/7rgsort/{sample}.part_042.clean.bam",
        p43="{path}preprocessing/7rgsort/{sample}.part_043.clean.bam",
        p44="{path}preprocessing/7rgsort/{sample}.part_044.clean.bam",
        p45="{path}preprocessing/7rgsort/{sample}.part_045.clean.bam",
        p46="{path}preprocessing/7rgsort/{sample}.part_046.clean.bam",
        p47="{path}preprocessing/7rgsort/{sample}.part_047.clean.bam",
        p48="{path}preprocessing/7rgsort/{sample}.part_048.clean.bam",
        p49="{path}preprocessing/7rgsort/{sample}.part_049.clean.bam",
        p50="{path}preprocessing/7rgsort/{sample}.part_050.clean.bam",
        p51="{path}preprocessing/7rgsort/{sample}.part_051.clean.bam",
        p52="{path}preprocessing/7rgsort/{sample}.part_052.clean.bam",
        p53="{path}preprocessing/7rgsort/{sample}.part_053.clean.bam",
        p54="{path}preprocessing/7rgsort/{sample}.part_054.clean.bam",
        p55="{path}preprocessing/7rgsort/{sample}.part_055.clean.bam",
        p56="{path}preprocessing/7rgsort/{sample}.part_056.clean.bam",
        p57="{path}preprocessing/7rgsort/{sample}.part_057.clean.bam",
        p58="{path}preprocessing/7rgsort/{sample}.part_058.clean.bam",
        p59="{path}preprocessing/7rgsort/{sample}.part_059.clean.bam",
        p60="{path}preprocessing/7rgsort/{sample}.part_060.clean.bam",
        p61="{path}preprocessing/7rgsort/{sample}.part_061.clean.bam",
        p62="{path}preprocessing/7rgsort/{sample}.part_062.clean.bam",
        p63="{path}preprocessing/7rgsort/{sample}.part_063.clean.bam",
        p64="{path}preprocessing/7rgsort/{sample}.part_064.clean.bam",
        p65="{path}preprocessing/7rgsort/{sample}.part_065.clean.bam",
        p66="{path}preprocessing/7rgsort/{sample}.part_066.clean.bam",
        p67="{path}preprocessing/7rgsort/{sample}.part_067.clean.bam",
        p68="{path}preprocessing/7rgsort/{sample}.part_068.clean.bam",
        p69="{path}preprocessing/7rgsort/{sample}.part_069.clean.bam",
        p70="{path}preprocessing/7rgsort/{sample}.part_070.clean.bam",
        p71="{path}preprocessing/7rgsort/{sample}.part_071.clean.bam",
        p72="{path}preprocessing/7rgsort/{sample}.part_072.clean.bam",
        p73="{path}preprocessing/7rgsort/{sample}.part_073.clean.bam",
        p74="{path}preprocessing/7rgsort/{sample}.part_074.clean.bam",
        p75="{path}preprocessing/7rgsort/{sample}.part_075.clean.bam",
        p76="{path}preprocessing/7rgsort/{sample}.part_076.clean.bam",
        p77="{path}preprocessing/7rgsort/{sample}.part_077.clean.bam",
        p78="{path}preprocessing/7rgsort/{sample}.part_078.clean.bam",
        p79="{path}preprocessing/7rgsort/{sample}.part_079.clean.bam",
        p80="{path}preprocessing/7rgsort/{sample}.part_080.clean.bam",
        p81="{path}preprocessing/7rgsort/{sample}.part_081.clean.bam",
        p82="{path}preprocessing/7rgsort/{sample}.part_082.clean.bam",
        p83="{path}preprocessing/7rgsort/{sample}.part_083.clean.bam",
        p84="{path}preprocessing/7rgsort/{sample}.part_084.clean.bam",
        p85="{path}preprocessing/7rgsort/{sample}.part_085.clean.bam",
        p86="{path}preprocessing/7rgsort/{sample}.part_086.clean.bam",
        p87="{path}preprocessing/7rgsort/{sample}.part_087.clean.bam",
        p88="{path}preprocessing/7rgsort/{sample}.part_088.clean.bam",
        p89="{path}preprocessing/7rgsort/{sample}.part_089.clean.bam",
        p90="{path}preprocessing/7rgsort/{sample}.part_090.clean.bam",
        p91="{path}preprocessing/7rgsort/{sample}.part_091.clean.bam",
        p92="{path}preprocessing/7rgsort/{sample}.part_092.clean.bam",
        p93="{path}preprocessing/7rgsort/{sample}.part_093.clean.bam",
        p94="{path}preprocessing/7rgsort/{sample}.part_094.clean.bam",
        p95="{path}preprocessing/7rgsort/{sample}.part_095.clean.bam",
        p96="{path}preprocessing/7rgsort/{sample}.part_096.clean.bam",
        p97="{path}preprocessing/7rgsort/{sample}.part_097.clean.bam",
        p98="{path}preprocessing/7rgsort/{sample}.part_098.clean.bam",
        p99="{path}preprocessing/7rgsort/{sample}.part_099.clean.bam",
        p100="{path}preprocessing/7rgsort/{sample}.part_100.clean.bam"
    output:
        "{path}preprocessing/8merged/{sample}.lanemerge.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
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

## Purge PCR duplicate reads
rule STEP11_purgeduplicates:
    input:
        "{path}preprocessing/8merged/{sample}.lanemerge.bam"
    output:
        "{path}preprocessing/9dedup/{sample}.dp.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
    
## Filter reads for only uniquely mapping
rule STEP12_mapqfilter:
    input:
        "{path}preprocessing/9dedup/{sample}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}.u.bam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        1
    resources:
        run_time=lambda params, attempt: attempt * 4
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
    
## Move bam to parent directory and rename
rule STEP13_move_bam:
    input:
        "{path}preprocessing/10unique/{sample}.u.bam"
    output:
        "{path}bam/{sample}.bam"
    shell:
        """
        mv {wildcards.path}preprocessing/10unique/*.bam {wildcards.path}bam/{wildcards.sample}.bam
        touch {output}
        """

## Build the .bai index for the processed bam file
rule STEP14_build_bai_index:
    input:
        "{path}bam/{sample}.bam"
    output:
        "{path}bam/{sample}.bam.bai"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

########################################################################################################################################
#### FOOTPRINTING - SAMPLE MERGED PEAKS ##############################################################################################
########################################################################################################################################

#### 
rule FOOTPRINTING_sample_merged_generate_binding_sites:
    input:
        "{path}peaks/sample_merged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_merged/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_binding_sites.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "snakeResources/scripts/generateBindingSites.R"

#### Generate the insertion matrix for sample specific peaks ####
rule FOOTPRINTING_sample_merged_generate_insertion_matrix:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "snakeResources/scripts/atacFunctions.R",
        "{path}footprints/sample_merged/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_binding_sites.RData"
    output:
        "{path}footprints/sample_merged/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_insertion_matrix.complete"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/generateInsertionMatrix.R"

#### Raw footprinting analysis, sample-specific peaks ####
rule FOOTPRINTING_sample_merged_raw_analysis:
    input:
        "{path}footprints/sample_merged/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_insertion_matrix.complete",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_merged/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_raw_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/analyzeRawFootprint.R"

rule EXAPNDER_footprinting_analysis_sample_merged:
    input:
        expand("{{path}}footprints/sample_merged/raw/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.sample_merged_raw_footprint.RData", genename=config["geneNames"])
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_merged_footprint_analysis.complete"
    shell:
        "touch {output}"

rule FOOTPRINTING_sample_merged_aggregate_footprint_data:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_merged_footprint_analysis.complete",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_merged/aggregated/{sample}.rep{repnum}.ref{refgenome}.sample_merged_aggregated_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/aggregateFootprintData.R"
