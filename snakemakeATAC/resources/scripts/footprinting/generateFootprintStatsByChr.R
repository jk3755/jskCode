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
  
  #### Report ####
  cat("Spooling footprint analysis", "\n")
  cat("Input path:", inputPath, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output filepath for footprint data:", dataOutPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  
  #### Sample name for the dataframe
  sampleID <- paste0(sampleName, ".", sampleRep)
  cat("Sample ID:", sampleID, "\n")
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", dataOutPath, "\n")
  if (file.exists(dataOutPath)){
    
    cat("Output file already exists. Skipping", "\n")
    
  } else {
    
    cat("Output file not found. Processing", "\n")
    
    #### Generate the input directory path ####
    inputDirectory <- str_replace(inputPath, "(?<=insertions/).*", paste0(geneName, "/"))
    cat("Input directory at path:", inputDirectory, "\n")
    
    #### Output directory
    outputDirectory <- str_replace(dataOutPath, "(?<=stats/).*", paste0(geneName, "/"))
    cat("Generating output directory at path:", outputDirectory, "\n")
    dir.create(path = outputDirectory, showWarnings = FALSE)
    
    #### Get the input file paths ####
    inputFileList <- list.files(path = inputDirectory, full.names = TRUE)
    numInputFiles <- length(inputFileList)
    cat("Identified", numInputFiles, "input files", "\n")
    
    #### If no input files, just make a dummy file ####
    if (numInputFiles == 0){
      
      cat("No input files. Making dummy output file", "\n")
      footprintSiteStatistics <- list()
      save(footprintSiteStatistics, file = dataOutPath)
      
    } else {
      
      #### Source functions ####
      cat("Loading functions from path:", functionSourcePath, "\n")
      source(functionSourcePath)
      
      #### Analyze the footprints ####
      for (a in 1:numInputFiles){
        
        ## Load the data
        cat("Loading data from path:", inputFileList[a], "\n")
        load(file = inputFileList[a])
        bindingSites <- insertionMatrixData$bindingSites
        insertionMatrix <- insertionMatrixData$insertionMatrix
        
        ## Get the scope for current file
        currentChr <- str_extract(inputFileList[a], "(?<=chr)([\\d]+|[\\w])")
        currentChr <- paste0("chr", currentChr)
        
        ##
        siteStatistics <- generateFootprintStats(insertionMatrix, bindingSites, sampleID, geneName)
        outputPath <- paste0(outputDirectory, currentChr, ".siteStats.", geneName, ".RData")
        save(siteStatistics, file = outputPath)
        
      }
 
      #### Touch the output file ####
      cat("Saving footprint data", "\n")
      file.create(dataOutPath)

    }
  }
}, finally = {
})