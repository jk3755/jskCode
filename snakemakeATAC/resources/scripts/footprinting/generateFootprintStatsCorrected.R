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
  
  ####
  seqbiasModelPath <- snakemake@input[[3]]
  seqbiasFastaPath <- snakemake@input[[4]]

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
        
        ##
        com <- paste0("SiteStatistics", a, " <- generateFootprintStatsSeqbiasCorrected(insertionMatrix, bindingSites, sampleID, seqbiasPath = seqbiasModelPath, fastaPath = seqbiasFastaPath)")
        eval(parse(text = com))
        
      }
      
      #### Merge all of the site statistics tables ####
      cat("Merging site statistics tables", "\n")
      currentTables <- as.character(paste0("SiteStatistics", 1:numInputFiles))
      currentTablesStr <- paste(currentTables, collapse = ",")
      com <- paste0("footprintSiteStatistics <- rbind(", currentTablesStr, ")")
      eval(parse(text = com))

      #### Save the footprint data ####
      cat("Saving footprint data", "\n")
      save(footprintSiteStatistics, file = dataOutPath)
    }
  }
}, finally = {
})