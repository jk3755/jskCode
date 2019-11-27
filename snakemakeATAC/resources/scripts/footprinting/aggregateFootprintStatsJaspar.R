## 
suppressMessages(library(stringr))

#### Set snakemake variables ####
cat("Setting snakemake variables", "\n")
inputPath <- snakemake@input[[1]]
functionSourcePath <- snakemake@input[[2]]
dataOutPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
geneName <- snakemake@wildcards[["ID"]]
refGenome <- snakemake@wildcards[["refgenome"]]
workingPath <- snakemake@wildcards[["path"]]

#### Report ####
cat("Aggregating footprint stats data", "\n")

#### Generate the input directory path ####
inputDirectory <- str_replace(inputPath, "(?<=stats/).*", paste0(geneName, "/"))
cat("Input directory at path:", inputDirectory, "\n")

#### Get the input file paths ####
inputFileList <- list.files(path = inputDirectory, full.names = TRUE)
numInputFiles <- length(inputFileList)
cat("Identified", numInputFiles, "input files", "\n")

##
if (numInputFiles == 0){
  
  cat("No input files. Making dummy output file", "\n")
  siteStatistics <- list()
  save(siteStatistics, file = dataOutPath)
  
} else {

  cat("Processing", "\n")
  ##
  load(inputFileList[1])
  aggregatedFootprintStats <- siteStatistics
  
  ####
  for (a in 2:numInputFiles){
    cat("File:", a, "\n")
    currentFile <- inputFileList[a]
    load(currentFile)
    aggregatedFootprintStats <- rbind(aggregatedFootprintStats, siteStatistics)
  }
  
  ####
  cat("Saving", "\n")
  save(aggregatedFootprintStats, file = dataOutPath)
  RDSout <- gsub("RData", "RDS", dataOutPath)
  saveRDS(aggregatedFootprintStats, file = RDSout)
}