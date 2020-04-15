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
  pwmScanScore <- snakemake@wildcards[["matchScore"]]
  
  #### Report ####
  cat("Generating binding sites:", "\n")
  cat("Peaks file:", peakPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("PWM matching score:", pwmScanScore, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output filepath for binding sites:", bindingSitesOutPath, "\n")
  
  #### Load libraries ####
  cat("Loading libraries", "\n")

  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  #### Retrieve the binding sites for the current gene ####
  cat("Retrieving binding sites", "\n")
  tempAllSites <- getAllBindingSites(geneName, pwmScanScore = pwmScanScore)
  
  #### Subset the binding sites with the peaks
  cat("Subsetting sites to peak regions", "\n")
  peaksGR <- importBED(peakPath)
  bindingSites <- subsetByOverlaps(tempAllSites, peaksGR)
  numSites <- length(bindingSites@ranges)
  cat("Identified", numSites, "total binding sites", "\n")
  
  #### Save the footprint data ####
  cat("Saving binding sites data", "\n")
  save(bindingSites, file = bindingSitesOutPath)
}, finally = {
})