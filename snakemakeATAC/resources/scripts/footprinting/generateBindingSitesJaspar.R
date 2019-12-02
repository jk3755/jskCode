#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  peakPath <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  jasparPath <- snakemake@input[[3]]
  bindingSitesOutPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneID <- snakemake@wildcards[["ID"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  pwmScanScore <- snakemake@wildcards[["matchScore"]]
  
  #### Report ####
  cat("Generating binding sites:", "\n")
  cat("Peaks file:", peakPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneID, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("PWM matching score:", pwmScanScore, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Filepath of JASPAR file:", jasparPath, "\n")
  cat("Output filepath for binding sites:", bindingSitesOutPath, "\n")
  
  ## Absolute filepath
  workingDir <- getwd()
  bindingSitesOutPath <- paste0(workingDir, "/", bindingSitesOutPath)
  cat("Output filepath for binding sites:", bindingSitesOutPath, "\n")
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", bindingSitesOutPath, "\n")
  if (file.exists(bindingSitesOutPath)){
    cat("Output file already exists. Skipping", "\n")
  } else {
    cat("Output file not found. Processing", "\n")
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    ##
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    suppressMessages(library(Biostrings))
    suppressMessages(library(MotifDb))
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(Rsamtools))
    suppressMessages(library(GenomicAlignments))
    suppressMessages(library(ChIPseeker))
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    suppressMessages(library(genomation))
    suppressMessages(library(TFBSTools))
    genome <- Hsapiens
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    ## Load JASPAR file
    load(file = jasparPath)
    
    ## Get binding sites and add score2
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
    
    ## Remove chrM entries
    cat("Removing mitochondrial binding sites, if present", "\n")
    indexChrM <- which(seqnames(allSites) == "chrM")
    if (length(indexChrM > 0)){
      allSites <- allSites[-indexChrM]
    }
    
    ## Perform final trim to standard only, reorder
    cat("Trimming to standard chromosomes only", "\n")
    allSites <- keepStandardChromosomes(allSites, pruning.mode = "coarse")
    allSites <- trim(allSites)
    allSites <- sortSeqlevels(allSites)
    allSites <- sort(allSites)
    
    ##
    peaksGR <- importBED(peakPath)
    bindingSites <- subsetByOverlaps(allSites, peaksGR)
    
    #### Save the footprint data ####
    cat("Saving binding sites data", "\n")
    save(bindingSites, file = bindingSitesOutPath)
  }
}, finally = {
})