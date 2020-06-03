#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  ## Prevent R from converting strings to factors (as is done in data frames)
  options(stringsAsFactors = FALSE)
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  functionSourcePath <- snakemake@input[[2]]
  dataOutPath <- snakemake@output[[1]]
  sampleDirectory <- snakemake@wildcards[["path"]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", dataOutPath, "\n")
  if (file.exists(dataOutPath)){
    cat("Output file already exists. Skipping", "\n")
  } else {
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    #### Input directory
    sampleString <- paste0(sampleName, ".rep", sampleRep, ".ref", refGenome, ".fp.RData")
    matchString <- paste0("aggregated/", sampleString)
    inputDirectory <- gsub(matchString, "stats/", dataOutPath)
    cat("Input directory at path:", inputDirectory, "\n")
    inputFileList <- list.files(path = inputDirectory, full.names = TRUE)
    numInputFiles <- length(inputFileList)
    cat("Identified", numInputFiles, "input files", "\n")
    
    #### Get the gene names ####
    geneNames <- c()
    
    for (a in 1:numInputFiles){
      currentFile <- inputFileList[a]
      geneName <- gsub(inputDirectory, "", currentFile)
      sampleString2 <- paste0(sampleName, ".rep", sampleRep, ".ref", refGenome, ".")
      geneName <- gsub(sampleString2, "", geneName)
      geneName <- gsub(".fp.RData", "", geneName)
      geneNames[a] <- geneName
    }
    
    ##
    footprintStatsList <- aggregateFootprintStats(inputFileList, geneNames)
    
    #### Save the aggregated footprint data ####
    cat("Saving aggregated footprint data", "\n")
    save(footprintStatsList, file = dataOutPath)
    
  }
}, finally = {
})
