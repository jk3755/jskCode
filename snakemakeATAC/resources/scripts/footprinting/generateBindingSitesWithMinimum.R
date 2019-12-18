#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  peakPath <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  bindingSitesOutPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  minimumSites <- snakemake@wildcards[["minsites"]]
  
  #### Report ####
  cat("Generating binding sites with minimum number of sites:", "\n")
  cat("Peaks file:", peakPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("PWM matching score:", pwmScanScore, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output filepath for binding sites:", bindingSitesOutPath, "\n")
  
  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  #### Retrieve the binding sites for the current gene ####
  cat("Retrieving binding sites", "\n")
  tempAllSites <- getAllBindingSitesWithMinimum(geneName, minimumSites, peakPath)
  numSites <- length(tempAllSites@ranges)
  cat("Identified", numSites, "total binding sites", "\n")
  
  #### Save the footprint data ####
  cat("Saving binding sites data", "\n")
  save(tempAllSites, file = bindingSitesOutPath)

}, finally = {
}) # end tryCatch