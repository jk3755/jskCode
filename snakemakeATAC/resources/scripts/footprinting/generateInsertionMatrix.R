#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  bamPath <- snakemake@input[[1]]
  baiPath <- snakemake@input[[2]]
  bindingSitesPath <- snakemake@input[[3]]
  functionSourcePath <- snakemake@input[[4]]
  dataOutPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  currentChunk <- as.numeric(snakemake@wildcards[["chunk"]])
  totalChunk <- as.numeric(snakemake@wildcards[["totalchunk"]])
  
  #### Report ####
  cat("Generating insertion matrix with the following parameters:", "\n")
  cat("Bam file:", bamPath, "\n")
  cat("Bai file:", baiPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("Current chunk:", currentChunk, "\n")
  cat("Total chunk:", totalChunk, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Filepath for loading binding sites:", bindingSitesPath, "\n")
  cat("Filepath for output data:", dataOutPath, "\n")
  
  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  #### Retrieve the binding sites for the current gene ####
  cat("Loading binding sites", "\n")
  load(bindingSitesPath)
  numSites <- length(bindingSites)
  
  #### No binding sites
  if (numSites == 0){
    
    cat("No binding sites identified, creating dummy file", "\n")
    insertionMatrixData <- "DUMMY"
    save(insertionMatrixData, file = dataOutPath)
    
  } else {
    
    #### Subset binding sites to current chunk ####
    cat("Subsetting binding sites", "\n")
    currentSites <- splitGRangesToBin(bindingSites, totalChunk, currentChunk)
    
    if (is.null(currentSites) | length(currentSites) == 0){
      cat("No binding sites in current chunk, outputting dummy file", "\n")
      insertionMatrixData <- "DUMMY"
      save(insertionMatrixData, file = dataOutPath)
      
    } else {
      
      ####
      numCurrentSites <- length(currentSites@ranges)
      cat("Identified", numCurrentSites, "in current chunk", "\n")
      maxWidth <- max(bindingSites@ranges@width)
      cat("Maximum width of binding sites:", maxWidth, "\n")
      
      ## Generate the matrix
      cat("Generating insertion matrix", "\n")
      insMatrix <- generateInsertionMatrix(bamPath, currentSites, maxWidth)
      
      ## Make the list
      cat("Storing data", "\n")
      insertionMatrixData <- list()
      insertionMatrixData$bindingSites <- currentSites
      insertionMatrixData$insertionMatrix <- insMatrix
      
      ## Save the file
      cat("Saving insertion matrix", "\n")
      save(insertionMatrixData, file = dataOutPath)
    }
  }
}, finally = {
})