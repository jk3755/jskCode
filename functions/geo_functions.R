
### Notes ####################################################################################################################################

## Retrieve data from article "Lineage-specific and single-cell chromatin accessibility charts human hematopoiesis and leukemia evolution"
## Corces et. al. Published: 15 August 2016, Nature
## https://www.nature.com/articles/ng.3646
## The data in this paper is available from the Gene Expression Omnibus (GEO)
## For info on how to query GEO data using R, see: http://genomicsclass.github.io/book/pages/GEOquery.html
## For info on how to query GEO, see: https://www.ncbi.nlm.nih.gov/geo/info/download.html


### Accession codes ##########################################################################################################################
## BioProject: PRJNA299657
## SRA: SRP065190
## GEO SuperSeries: GSE75384


### Code #####################################################################################################################################
## Install libraries, if necessary
#chooseCRANmirror()
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("GEOquery")

## Load libraries
library(GEOquery)
library(GenomicRanges)

## Set accession code
accID <- "GSE75384"

## Query the data
allData <- getGEO(accID, GSEMatrix = TRUE)

## Split up the data
rnaData <- allData[["GSE75384-GPL11154_series_matrix.txt.gz"]]
rnaData2 <-allData[["GSE75384-GPL18573_series_matrix.txt.gz"]]
atacData <- allData[["GSE75384-GPL16791_series_matrix.txt.gz"]]


#########################################################################################################
#### Get the ATAC-seq peaks counts for all samples
# file origin: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74912
# ftp link: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz
# http link: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74310&format=file&file=GSE74310%5FscATACseq%5FAll%5FCounts%2Etxt%2Egz
# ftpATAC <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz"
# httpATAC <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74310&format=file&file=GSE74310%5FscATACseq%5FAll%5FCounts%2Etxt%2Egz"
# url <- getURL(httpATAC)

##
atacPeaksPath <- "C:\\Users\\jsk33\\OneDrive\\DATA_COSMA\\corces2016data\\GSE74912_ATACseq_All_Counts.txt"
atacPeaksCounts <- read.table(atacPeaksPath, header = TRUE)
atacSampleNames <- names(atacPeaksCounts)
atacSampleNames <- atacSampleNames[4:135]

## Write to file
atacPeaksOutpath <- "C:\\Users\\jsk33\\OneDrive\\DATA_COSMA\\corces2016data\\GSE74912_ATACseq_All_Counts.RData"
save(atacPeaksCounts, file = atacPeaksOutpath)
load(atacPeaksOutpath)

## Write peaks to BED file
sourcePath <- "C:\\Users\\jsk33\\OneDrive\\3.Git\\jskCode\\snakemakeATAC\\resources\\functions\\atacFunctions.R"
source(sourcePath)
atacPeaksDF <- data.frame(chr = atacPeaksCounts[,1], start = atacPeaksCounts[,2], end = atacPeaksCounts[,3])
atacPeaksGR <- makeGRangesFromDataFrame(atacPeaksDF)
peaksOutpath <- "C:\\Users\\jsk33\\OneDrive\\DATA_COSMA\\corces2016data\\GSE74912_ATACseq_All_Counts.narrowPeak"
writeGRangesToBED(atacPeaksGR, peaksOutpath)
