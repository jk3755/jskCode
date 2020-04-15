---
  title: "cosmaMergedAnalysis"
author: "Jordan S. Kesner"
date: "1/14/2020"
output: pdf_document
editor_options: 
  chunk_output_type: console
---
  
  Setup
```{r}

#### Load ####
suppressPackageStartupMessages({
  library(enrichR)
  library(genomation)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(enrichR)
  library(ChIPseeker)
  library(org.Hs.eg.db)
  library(biomaRt)
  #library(xlsx)
  library(BiocParallel)
  library(Rsamtools)
  library(ShortRead)
  library(RIdeogram)
  library(chromVAR)
  library(motifmatchr)
  library(Matrix)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(motifmatchr)
  library(chromVARmotifs)
  library(cowplot)
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
})

#### Options ####
options(scipen = 999)
seed = 666

#### Parallelilzation ####
register(SnowParam(workers = 12, type = "SOCK")) # for windows
#register(MulticoreParam(4)) # Use 8 cores # for linux

#### Set onedrive path ####
path.onedrive <- "C:\\Users\\jsk33\\OneDrive\\" ## home computer
#path.onedrive <- "/home/local/ARCS/jk3755/OneDrive/" ## lab linux
#path.onedrive <- "C:\\Users\\Jordan\\OneDrive\\" ## lab laptop

#### Source functions ####
path.granges_functions <- file.path(path.onedrive, "3.Git/jskCode/functions/granges_functions.R")
source(path.granges_functions)

```

Set filepaths to data files
```{r}

#### Set BAM paths for our samples ####
path.bam.baz2b_A <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/bam/DonorA_Baz2B.rep1.refhg38.bam")
path.bam.luciferase_A <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/bam/DonorA_Luf.rep1.refhg38.bam")
path.bam.progenitor_A <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/bam/DonorA_Progenitor.rep1.refhg38.bam")
path.bam.baz2b_B <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/bam/DonorB_Baz2B.rep1.refhg38.bam")
path.bam.luciferase_B <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/bam/DonorB_Luf.rep1.refhg38.bam")
path.bam.progenitor_B <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/bam/DonorB_Progenitor.rep1.refhg38.bam")

#### Paths for preliminary chromVar analysis with 12 Corces2016 samples ####
#### 4983 ####
path.bam.4983.HSC <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920466.bam")
path.bam.4983.MPP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920467.bam")
path.bam.4983.LMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920468.bam")
path.bam.4983.CMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920469.bam")
path.bam.4983.GMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920471.bam")
path.bam.4983.MEP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920473.bam")

#### 6792 ####
path.bam.6792.HSC <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920477.bam")
path.bam.6792.MPP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920479.bam")
path.bam.6792.LMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920480.bam")
path.bam.6792.CMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920481.bam")
path.bam.6792.GMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920483.bam")
path.bam.6792.MEP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/bam/SRR2920485.bam")

####
path.peaks.4983.HSC <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920466_globalnorm_p001_peaks.narrowPeak")
path.peaks.4983.MPP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920467_globalnorm_p001_peaks.narrowPeak")
path.peaks.4983.LMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920468_globalnorm_p001_peaks.narrowPeak")
path.peaks.4983.CMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920469_globalnorm_p001_peaks.narrowPeak")
path.peaks.4983.GMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920471_globalnorm_p001_peaks.narrowPeak")
path.peaks.4983.MEP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920473_globalnorm_p001_peaks.narrowPeak")

####
path.peaks.6792.HSC <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920477_globalnorm_p001_peaks.narrowPeak")
path.peaks.6792.MPP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920479_globalnorm_p001_peaks.narrowPeak")
path.peaks.6792.LMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920480_globalnorm_p001_peaks.narrowPeak")
path.peaks.6792.CMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920481_globalnorm_p001_peaks.narrowPeak")
path.peaks.6792.GMP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920483_globalnorm_p001_peaks.narrowPeak")
path.peaks.6792.MEP <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/peaks/SRR2920485_globalnorm_p001_peaks.narrowPeak")

```

Make merged peak set from Corces analysis
```{r}

##
peaks.4983.HSC <- importNarrowpeakAsGRangesStandardized(path.peaks.4983.HSC)
peaks.4983.MPP <- importNarrowpeakAsGRangesStandardized(path.peaks.4983.MPP)
peaks.4983.LMP <- importNarrowpeakAsGRangesStandardized(path.peaks.4983.LMP)
peaks.4983.CMP <- importNarrowpeakAsGRangesStandardized(path.peaks.4983.CMP)
peaks.4983.GMP <- importNarrowpeakAsGRangesStandardized(path.peaks.4983.GMP)
peaks.4983.MEP <- importNarrowpeakAsGRangesStandardized(path.peaks.4983.MEP)

##
peaks.6792.HSC <- importNarrowpeakAsGRangesStandardized(path.peaks.6792.HSC)
peaks.6792.MPP <- importNarrowpeakAsGRangesStandardized(path.peaks.6792.MPP)
peaks.6792.LMP <- importNarrowpeakAsGRangesStandardized(path.peaks.6792.LMP)
peaks.6792.CMP <- importNarrowpeakAsGRangesStandardized(path.peaks.6792.CMP)
peaks.6792.GMP <- importNarrowpeakAsGRangesStandardized(path.peaks.6792.GMP)
peaks.6792.MEP <- importNarrowpeakAsGRangesStandardized(path.peaks.6792.MEP)

#### Put them in a list for merging ####
list.peaks <- list(peaks.4983.HSC, peaks.4983.MPP, peaks.4983.LMP, peaks.4983.CMP, peaks.4983.GMP, peaks.4983.MEP,
                   peaks.6792.HSC, peaks.6792.MPP, peaks.6792.LMP, peaks.6792.CMP, peaks.6792.GMP, peaks.6792.MEP)

#### Merge peaks ####
message(">>> Merging peaks <<<")
merged_peaks <- unique(getGRangesUnion(list.peaks))

#### Clean the peaks ####
message(">>> Cleaning merged peaks <<<")
merged_peaks <- cleanGRanges(merged_peaks)

```

Merged peaks from our samples
```{r}

####
path.peaks.BAZ2B_A <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/peaks/BAZ2B_A_globalnorm_p001_peaks.narrowPeak")
path.peaks.LUF_A <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/peaks/LUF_A_globalnorm_p001_peaks.narrowPeak")
path.peaks.PROG_A <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/peaks/PROG_A_globalnorm_p001_peaks.narrowPeak")
path.peaks.BAZ2B_B <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/peaks/BAZ2B_B_globalnorm_p001_peaks.narrowPeak")
path.peaks.LUF_B <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/peaks/LUF_B_globalnorm_p001_peaks.narrowPeak")
path.peaks.PROG_B <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/peaks/PROG_B_globalnorm_p001_peaks.narrowPeak")

####
peaks.BAZ2B_A <- importNarrowpeakAsGRangesStandardized(path.peaks.BAZ2B_A)
peaks.LUF_A <- importNarrowpeakAsGRangesStandardized(path.peaks.LUF_A)
peaks.PROG_A <- importNarrowpeakAsGRangesStandardized(path.peaks.PROG_A)
peaks.BAZ2B_B <- importNarrowpeakAsGRangesStandardized(path.peaks.BAZ2B_B)
peaks.LUF_B <- importNarrowpeakAsGRangesStandardized(path.peaks.LUF_B)
peaks.PROG_B <- importNarrowpeakAsGRangesStandardized(path.peaks.PROG_B)

#### Put them in a list for merging ####
list.peaks <- list(peaks.BAZ2B_A, peaks.LUF_A, peaks.PROG_A,
                   peaks.BAZ2B_B, peaks.LUF_B, peaks.PROG_B)

#### Merge peaks ####
message(">>> Merging peaks <<<")
merged_peaks <- unique(getGRangesUnion(list.peaks))

#### Clean the peaks ####
message(">>> Cleaning merged peaks <<<")
merged_peaks <- cleanGRanges(merged_peaks)

```

```{r}

interesting_peaks <- readRDS(file = paste0(path.onedrive, "4.Analysis\\ATAC_COSMA\\data\\summits\\interesting_windows.RDS"))
merged_peaks <- cleanGRanges(interesting_peaks)

interesting_peaks <- readRDS(file = paste0(path.onedrive, "4.Analysis\\ATAC_COSMA\\data\\summits\\interesting_windows2.RDS"))
merged_peaks <- cleanGRanges(interesting_peaks)

```


chromVar analysis
```{r}

#### Set the peaks
input_peaks <- interesting_peaks

#### Output directory ####
dir.chromVar <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/chromVar/D/")

#### BAM filepath list ####
list.paths.bam <- c(path.bam.baz2b_A, path.bam.luciferase_A, path.bam.progenitor_A,
                    path.bam.baz2b_B, path.bam.luciferase_B, path.bam.progenitor_B,
                    path.bam.4983.HSC, path.bam.6792.HSC,
                    path.bam.4983.MPP, path.bam.6792.MPP,
                    path.bam.4983.LMP, path.bam.6792.LMP,
                    path.bam.4983.CMP, path.bam.6792.CMP,
                    path.bam.4983.GMP, path.bam.6792.GMP,
                    path.bam.4983.MEP, path.bam.6792.MEP)


#### Sample IDs ####
chromvar.sampleIDs <- c("BAZ2B", "LUF", "PROG",
                        "BAZ2B", "LUF", "PROG",
                        "HSC", "HSC",
                        "MPP", "MPP",
                        "LMP", "LMP",
                        "CMP", "CMP",
                        "GMP", "GMP",
                        "MEP", "MEP")

#### Column data ####
colData = DataFrame(condition = as.factor(c("BAZ2B", "LUF", "PROG",
                                            "BAZ2B", "LUF", "PROG",
                                            "HSC", "HSC",
                                            "MPP", "MPP",
                                            "LMP", "LMP",
                                            "CMP", "CMP",
                                            "GMP", "GMP",
                                            "MEP", "MEP")))

#### Sample types ####
sampleTypes <- as.factor(c("BAZ2B", "LUF", "PROG",
                           "BAZ2B", "LUF", "PROG",
                           "HSC", "HSC",
                           "MPP", "MPP",
                           "LMP", "LMP",
                           "CMP", "CMP",
                           "GMP", "GMP",
                           "MEP", "MEP"))

#### Detect # of samples ####
numSamples <- length(list.paths.bam)
message(">>> Number of samples detected: " , numSamples)

#### Generate the fragment counts ####
message(">>> Generating fragment counts <<<")
fragment_counts <- getCounts(list.paths.bam,
                             input_peaks,
                             paired = TRUE,
                             format = "bam")

#### Save the fragment counts file ####
saveRDS(fragment_counts, file = paste0(dir.chromVar, "/fragment_counts_D.RDS"))

#### Begin the analysis ####
message(">>> Analyzing data <<<")
message(">>> Setting seed to:", seed)
set.seed(seed)
message(">>> Total peaks: ", nrow(fragment_counts))

#### Add experimental condition ####
message(">>> Setting experimental condition metadata <<<")
fragment_counts@colData$condition <- colData

#### Adjust for GC bias ####
message(">>> Adjusting for GC Bias <<<")
fragment_counts <- addGCBias(fragment_counts,
                             genome = BSgenome.Hsapiens.UCSC.hg38)

#### Filter peaks ####
message(">>> Filtering peaks <<<")
fragment_counts <- sort(fragment_counts)
filtered_counts <- filterPeaks(fragment_counts, non_overlapping = FALSE)
message(">>> Total peaks after filtering: ", nrow(filtered_counts))

#### Load PWM motifs ####
message(">>> Loading PWM Motifs <<<")
data("human_pwms_v2",verbose = FALSE)
message(">>> Using " , length(human_pwms_v2), " motifs")

#### Match the motifs to the counts ####
message(">>> Matching the motifs <<<")
motif_ix <- matchMotifs(pwms = human_pwms_v2 ,
                        subject = filtered_counts, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

#### Compute the background peaks ####
message(">>> Computing background peak expectation <<<")
filtered_counts@colData$condition <- sampleTypes
expected <- computeExpectations(filtered_counts,
                                norm = TRUE,
                                group = colData(filtered_counts)$condition)
bg <- getBackgroundPeaks(object = filtered_counts)

#### Compute the deviation ####
message(">>> Computing deviation <<<")
dev <- computeDeviations(object = filtered_counts,
                         annotations = motif_ix,
                         expectation = expected,
                         background_peaks = bg)

#### Compute variability ####
message(">>> Computing variability <<<<")
variability <- computeVariability(dev)

####
makeTSNEvariablityPlot(variability, dev, 0.1, 5, "GATA2", file.path(dir.chromVar, "variability_plot.pdf"))


#################

#### Variability plot
message(">>> Variability Plot")
{
  my_threshold <- 7.5
  variability_table <- as_tibble(variability)
  variability_table$my_score <- variability_table$variability
  
  variability_table <- variability_table %>% arrange(desc(my_score))
  variability_table$name <- factor( variability_table$name , levels = variability_table$name , ordered = TRUE )
  variability_table$motif_order <- 1:nrow(variability_table)
  
  my_plot <- ggplot( data = variability_table , aes( x = motif_order , y = my_score ) ) +
    
    geom_point( data = variability_table %>% filter( my_score < my_threshold ) , alpha = 0.5 , stroke = 1 , size = 1 , color = "gray" , shape = 20 ) +
    geom_point( data = variability_table %>% filter( my_score >= my_threshold ) , alpha = 0.5 , stroke = 1 , size = 2 , color = "orange" , shape = 20 ) +
    
    geom_text_repel( data = variability_table %>% filter( my_score >= my_threshold ) , aes( label = name ), 
                     size = 3 , 
                     # point.padding = 0.2 , 
                     segment.alpha = 0.5 , segment.colour = "gray" ) +
    
    xlab( paste0( "Scored Motifs") ) +
    ylab( paste0( "log10 Variability") ) +
    ggtitle( paste0( "Variability Plot") ) +
    ylim(0, 30)
  theme_minimal() +
    theme( 
      plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.x = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.y = element_text(size=15,face="bold") ) # family="Open Sans") )
  
  print(my_plot)
  pdf(file.path(dir.chromVar,"variability-plot_withCorces2016samples.pdf") )
  print(my_plot)
  dev.off()
}







```


```{r}



#### Load the Corces peaks ####
#path.peaks.corces2016 <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/counts/GSE74912_ATACseq_All_Counts.narrowPeak")
path.peakCounts.corces2016 <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/data/corces2016/counts/GSE74912_ATACseq_All_Counts.RData")
load(path.peakCounts.corces2016)
df.corces2016.peaks <- data.frame(chr = atacPeaksCounts$Chr, start = atacPeaksCounts$Start, end = atacPeaksCounts$End)
peaks.corces2016 <- makeGRangesFromDataFrame(df.corces2016.peaks)

#### Import peak files ####
#peaks.BAZ2B_A <- importNarrowpeakAsGRangesStandardized(path.peak.baz2b_A)
#peaks.LUF_A <- importNarrowpeakAsGRangesStandardized(path.peak.luciferase_A)
#peaks.PROG_A <- importNarrowpeakAsGRangesStandardized(path.peak.progenitor_A)
#peaks.BAZ2B_B <- importNarrowpeakAsGRangesStandardized(path.peak.baz2b_B)
#peaks.LUF_B <- importNarrowpeakAsGRangesStandardized(path.peak.luciferase_B)
#peaks.PROG_B <- importNarrowpeakAsGRangesStandardized(path.peak.progenitor_B)

#### Put them in a list for merging ####
#list.peaks <- list(peaks.BAZ2B_A, peaks.LUF_A, peaks.PROG_A,
#                   peaks.BAZ2B_B, peaks.LUF_B, peaks.PROG_B)

#### Merge peaks ####
#message(">>> Merging peaks <<<")
#merged_peaks <- unique(getGRangesUnion(list.peaks))

#### Clean the peaks ####
#message(">>> Cleaning merged peaks <<<")
#merged_peaks <- cleanGRanges(merged_peaks)

#### Plot variability ####
message(">>> Plotting variability <<<<")
{
  p0 <- plotVariability(variability, use_plotly = FALSE)
  sample_cor <- getSampleCorrelation(dev)
  set.seed(seed)
  tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 2)
  tsne_plots <- plotDeviationsTsne(dev, tsne_results, 
                                   annotation_name = "SP1", 
                                   sample_column = "condition", 
                                   shiny = FALSE)
  
  p1 <- tsne_plots[[1]]
  p2 <- tsne_plots[[2]]
  plot_bottom_row <- plot_grid(p1, p2,ncol=2)
  p <- plot_grid(p0, plot_bottom_row, 
                 labels = c('A', 'B'),
                 nrow = 2, ncol = 1)
  
  pdf(file.path(dir.chromVar, "deviations-tsne.pdf") )
  print(p)
  dev.off()	
}
```


chromVar analysis without prog samples
```{r}

#### Set the peaks
input_peaks <- readRDS(paste0(path.onedrive, "4.Analysis\\ATAC_COSMA\\data\\summits\\distal_interesting_ranges.RDS"))

#### Output directory ####
dir.chromVar <- file.path(path.onedrive, "4.Analysis/ATAC_COSMA/chromVar/F/")

#### BAM filepath list ####
list.paths.bam <- c(path.bam.baz2b_A, path.bam.luciferase_A,
                    path.bam.baz2b_B, path.bam.luciferase_B,
                    path.bam.4983.HSC, path.bam.6792.HSC,
                    path.bam.4983.MPP, path.bam.6792.MPP,
                    path.bam.4983.LMP, path.bam.6792.LMP,
                    path.bam.4983.CMP, path.bam.6792.CMP,
                    path.bam.4983.GMP, path.bam.6792.GMP,
                    path.bam.4983.MEP, path.bam.6792.MEP)


#### Sample IDs ####
chromvar.sampleIDs <- c("BAZ2B", "LUF",
                        "BAZ2B", "LUF",
                        "HSC", "HSC",
                        "MPP", "MPP",
                        "LMP", "LMP",
                        "CMP", "CMP",
                        "GMP", "GMP",
                        "MEP", "MEP")

#### Column data ####
colData = DataFrame(condition = as.factor(c("BAZ2B", "LUF",
                                            "BAZ2B", "LUF",
                                            "HSC", "HSC",
                                            "MPP", "MPP",
                                            "LMP", "LMP",
                                            "CMP", "CMP",
                                            "GMP", "GMP",
                                            "MEP", "MEP")))

#### Sample types ####
sampleTypes <- as.factor(c("BAZ2B", "LUF",
                           "BAZ2B", "LUF",
                           "HSC", "HSC",
                           "MPP", "MPP",
                           "LMP", "LMP",
                           "CMP", "CMP",
                           "GMP", "GMP",
                           "MEP", "MEP"))

#### Detect # of samples ####
numSamples <- length(list.paths.bam)
message(">>> Number of samples detected: " , numSamples)

#### Generate the fragment counts ####
message(">>> Generating fragment counts <<<")
fragment_counts <- getCounts(list.paths.bam,
                             input_peaks,
                             paired = TRUE,
                             format = "bam")

#### Save the fragment counts file ####
saveRDS(fragment_counts, file = paste0(dir.chromVar, "/fragment_counts_F.RDS"))

#### Begin the analysis ####
message(">>> Analyzing data <<<")
message(">>> Setting seed to:", seed)
set.seed(seed)
message(">>> Total peaks: ", nrow(fragment_counts))

#### Add experimental condition ####
message(">>> Setting experimental condition metadata <<<")
fragment_counts@colData$condition <- colData

#### Adjust for GC bias ####
message(">>> Adjusting for GC Bias <<<")
fragment_counts <- addGCBias(fragment_counts,
                             genome = BSgenome.Hsapiens.UCSC.hg38)

#### Filter peaks ####
message(">>> Filtering peaks <<<")
fragment_counts <- sort(fragment_counts)
filtered_counts <- filterPeaks(fragment_counts, non_overlapping = FALSE)
message(">>> Total peaks after filtering: ", nrow(filtered_counts))

#### Load PWM motifs ####
message(">>> Loading PWM Motifs <<<")
data("human_pwms_v2",verbose = FALSE)
message(">>> Using " , length(human_pwms_v2), " motifs")

#### Match the motifs to the counts ####
message(">>> Matching the motifs <<<")
motif_ix <- matchMotifs(pwms = human_pwms_v2 ,
                        subject = filtered_counts, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

#### Compute the background peaks ####
message(">>> Computing background peak expectation <<<")
filtered_counts@colData$condition <- sampleTypes
expected <- computeExpectations(filtered_counts,
                                norm = TRUE,
                                group = colData(filtered_counts)$condition)
bg <- getBackgroundPeaks(object = filtered_counts)

#### Compute the deviation ####
message(">>> Computing deviation <<<")
dev <- computeDeviations(object = filtered_counts,
                         annotations = motif_ix,
                         expectation = expected,
                         background_peaks = bg)

#### Compute variability ####
message(">>> Computing variability <<<<")
variability <- computeVariability(dev)

####
makeTSNEvariablityPlot(variability, dev, 1, 4, "GATA2", file.path(dir.chromVar, "variability_plot.pdf"))

```




