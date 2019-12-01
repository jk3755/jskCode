#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
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
  currentChr <- snakemake@wildcards[["chr"]]
  
  #### Report ####
  cat("Spooling footprint analysis", "\n")
  cat("Input path:", inputPath, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output filepath for footprint data:", dataOutPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("Current chromosome:", currentChr, "\n")
  
  #### Sample name for the dataframe
  sampleID <- paste0(sampleName, ".", sampleRep)
  cat("Sample ID:", sampleID, "\n")

  #### Load the input file ####
  cat("Loading input file", "\n")
  load(inputPath)
  
  #### Check if dummy file ####
  if (insertionMatrixData == "DUMMY"){
    
    cat("Dummy file detected. Generating new dummy file", "\n")
    footprintSiteStatistics <- "DUMMY"
    save(footprintSiteStatistics, file = dataOutPath)
    
  } else {
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    #### Analyze the footprints ####
    cat("Analyzing footprints", "\n")
    bindingSites <- insertionMatrixData$bindingSites
    insertionMatrix <- insertionMatrixData$insertionMatrix
    siteStatistics <- generateFootprintStats(insertionMatrix, bindingSites, sampleID, geneName)
    save(siteStatistics, file = dataOutPath)
      
    }

}, finally = {
})
