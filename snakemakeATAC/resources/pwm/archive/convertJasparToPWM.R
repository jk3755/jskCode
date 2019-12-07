
## raw data
rawDir <- "C:\\Users\\jsk33\\OneDrive\\3.Git\\jskCode\\snakemakeATAC\\resources\\jaspar\\raw"

## output dir
outputDir <- "C:\\Users\\jsk33\\OneDrive\\3.Git\\jskCode\\snakemakeATAC\\resources\\jaspar\\pwm"
  
## source functions
source("C:\\Users\\jsk33\\OneDrive\\3.Git\\jskCode\\snakemakeATAC\\resources\\functions\\atacFunctions.R")

## input files
inputFileList <- list.files(path = rawDir, full.names = TRUE)
numInputFiles <- length(inputFileList)

## Load required libraries
suppressMessages(library(TFBSTools))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(GenomicRanges))
genome <- Hsapiens


### lets try this a different way
for (a in 1:numInputFiles){
  
  ## Load JASPAR file
  matrix <- readJASPARMatrix(inputFileList[a], matrixClass = "PWM")
  geneName <- paste0("KARTHIK", a)
  PWM <- matrix@listData[[1]]@profileMatrix
  outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
  save(PWM, file = outPath)
  
}







############################
for (a in 46:numInputFiles){
  
  ## Load JASPAR file
  matrix <- readJASPARMatrix(inputFileList[a], matrixClass = "PWM")
  geneName <- matrix@listData[[1]]@name
  PWM <- matrix@listData[[1]]@profileMatrix
  outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
  save(PWM, file = outPath)
  
}

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[5], matrixClass = "PWM")
geneName <- "DDIT3"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[16], matrixClass = "PWM")
geneName <- "TAL1"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[21], matrixClass = "PWM")
geneName <- "NFIC"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[23], matrixClass = "PWM")
geneName <- "GATA1"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[38], matrixClass = "PWM")
geneName <- "NR1H3"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[43], matrixClass = "PWM")
geneName <- "STAT1"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)

## some need manual curation
matrix <- readJASPARMatrix(inputFileList[45], matrixClass = "PWM")
geneName <- "BACH1"
PWM <- matrix@listData[[1]]@profileMatrix
outPath <- paste0(outputDir, "\\", geneName, ".pwm.RData")
save(PWM, file = outPath)