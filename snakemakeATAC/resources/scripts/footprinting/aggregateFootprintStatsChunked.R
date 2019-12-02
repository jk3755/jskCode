## 
suppressMessages(library(stringr))

#### Set snakemake variables ####
cat("Setting snakemake variables", "\n")
inputPath <- snakemake@input[[1]]
functionSourcePath <- snakemake@input[[2]]
dataOutPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
geneName <- snakemake@wildcards[["gene"]]
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

#### Initialize the empty data frame
cat("Initializing data frame", "\n")
aggregatedFootprintStats <- data.frame(sampleID = as.character(c()),
                                       geneName = as.character(c()),
                                       chr = as.character(c()),
                                       start = as.numeric(c()),
                                       width = as.numeric(c()),
                                       score = as.numeric(c()),
                                       score2 = as.numeric(c()),
                                       totalSignal = as.numeric(c()),
                                       averageTotalSignal = as.numeric(c()),
                                       totalMotifSignal = as.numeric(c()),
                                       averageMotifSignal = as.numeric(c()),
                                       totalFlankSignal = as.numeric(c()),
                                       averageFlankSignal = as.numeric(c()),
                                       totalBackgroundSignal = as.numeric(c()),
                                       averageBackgroundSignal = as.numeric(c()),
                                       flankAccessibility = as.numeric(c()),
                                       footprintDepth = as.numeric(c()),
                                       log2FlankAccessibility = as.numeric(c()),
                                       log2FootprintDepth = as.numeric(c()),
                                       binding = as.numeric(c()),
                                       pvalue = as.numeric(c()),
                                       zscore = as.numeric(c()),
                                       annotation = as.character(c()),
                                       symbol = as.character(c()),
                                       stringsAsFactors = FALSE) 

#### Load the RDS files and bind to dataframe ####
cat("Loading and appending data", "\n")
for (a in 1:numInputFiles){
  
  load(inputFileList[a])
  
  if (siteStatistics == "DUMMY"){
  } else {
    aggregatedFootprintStats <- rbind(aggregatedFootprintStats, siteStatistics)
  }
}

####
cat("Saving aggregated data", "\n")
save(aggregatedFootprintStats, file = dataOutPath)
RDSout <- gsub("RData", "RDS", dataOutPath)
saveRDS(aggregatedFootprintStats, file = RDSout)
