

#### Load averaged bigwig from RDS file by chromosome ####
loadBigwigRDSbyChr <- function(filepath, chromosome){
  
  ## Load the whole file
  temp <- readRDS(filepath)
  
  ## Subset to desired chr
  subtemp <- temp[seqnames(temp) == chromosome]
  
  ## cleanup
  rm(temp)
  gc()
  
  ## return
  return(subtemp)
  
}

#### Build a list object with many filepaths
constructFilepathList <- function(sampleNames, fileIdentifiers, filePrefix, fileSuffix, parentDir){
  
  ## Detect number of samples
  numSamples <- length(sampleNames)
  cat("Detected", numSamples, "input samples \n")
  
  ##
  cat("Parent directory:", parentDir, "\n")
  cat("Filename prefix:", filePrefix, "\n")
  cat("Filename suffix:", fileSuffix, "\n")
  
  ##
  tempList <- list()
  
  ##
  for(a in 1:numSamples){
    
    ##
    tempPath <- file.path(parentDir, paste0(filePrefix, fileIdentifiers[a], fileSuffix))
    
    ##
    com <- paste0("tempList$", sampleNames[a], " <- tempPath")
    eval(parse(text = com))
    
  }
  
  ##
  return(tempList)
  
}

#### 
loadObjectsFromFilepathList <- function(sampleNames, objectPrefix, filepathList, importFunction){
  
  ##
  numSamples <- length(sampleNames)
  
  ##
  for(a in 1:numSamples){
    
    ## filepath
    com <- paste0("tempPath <- filepathList$", sampleNames[a])
    eval(parse(text = com))
    
    ## object name
    tempName <- paste0(objectPrefix, sampleNames[a])
    
    ## import
    tempObject <- importFunction(tempPath)
    
    ## push object to global
    com <- paste0(tempName, "<<- tempObject")
    eval(parse(text = com))
    
  }
}

####
convertBED10colTo3col <- function(inputBED, outputPath){
  
  inBED <- importBED(inputBED)
  
  chrVec <- as.character(inBED@seqnames)
  startVec <- inBED@ranges@start
  widthVec <- inBED@ranges@width
  endVec <- startVec + widthVec
  
  outBED <- data.frame(chr = chrVec, start = startVec, end = endVec)
  outBED <- makeGRangesFromDataFrame(outBED)
  
  writeGRangesToBED(outBED, outputPath)
}

#### Import a BED file as a GRanges object ####
importBED <- function(bedPath){
  
  ## Report
  cat("Running importBED function", "\n")
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(genomation))
  
  ## 
  cat("Importing BED data from path:", bedPath, "\n")
  bedin <- readBed(bedPath, track.line = "AUTO", remove.unusual = TRUE, zero.based = TRUE)
  
  ## 
  cat("Cleaning GRanges object", "\n")
  bedin <- cleanGRanges(bedin)
  
  ## 
  cat("Returning data", "\n")
  return(bedin)
  
} # end importBED

#### Import a BED narrowPeak file output from MACS2 ####
importMACS2narrowPeak <- function(bedPath){
  
  ## Report
  cat("Running importMACS2narrowPeak function", "\n")
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(genomation))
  
  ## 
  cat("Importing BED data from path:", bedPath, "\n")
  bedin <- readBed(bedPath, track.line = "AUTO", remove.unusual = TRUE, zero.based = TRUE)
  
  ## 
  cat("Cleaning GRanges object", "\n")
  bedin <- cleanGRanges(bedin)
  
  #### Need to modify the metadata column names ####
  mscore <- mcols(bedin)$score
  mname <- mcols(bedin)$name
  msummitfc <- mcols(bedin)$thickStart
  msummitneglog10pvalue <- mcols(bedin)$thickEnd
  msummitneglog10qvalue <- mcols(bedin)$itemRgb
  msummitposition <- mcols(bedin)$blockCount
  ##
  df <- data.frame(score = mscore,
                   name = mname,
                   summitfc = msummitfc,
                   summitneglog10pvalue = msummitneglog10pvalue,
                   summitneglog10qvalue = msummitneglog10qvalue,
                   summitposition = msummitposition)
  
  ##
  mcols(bedin) <- NULL
  
  ##
  mcols(bedin, use.names = FALSE) <- df
  
  ##
  cat("Returning data", "\n")
  return(bedin)
  
} # end importMACS2narrowPeak

#### Generate BED interval for a particular chromosome of hg38
generateBEDinterval <- function(chrName = "chr1", outputPath){
  
  ## Report
  cat("Running generateBEDinterval function", "\n")
  cat("Selected chromosome:", chrName, "\n")
  cat("Output path for BED file:", outputPath, "\n")
  
  ## Load required libraries
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(genomation))
  
  ##
  cat("Loading BSgenome hg38 as GRanges", "\n")
  hg38GR <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
  
  ##
  cat("Subsetting to given chromosome", "\n")
  com <- paste0("refGR <- hg38GR[seqnames(hg38GR) == '", chrName, "']")
  eval(parse(text = com))
  refGR[2] <- refGR
  strand(refGR[1]) <- "+"
  strand(refGR[2]) <- "-"
  refGR@strand@values <- droplevels(refGR@strand@values)
  refGR <- keepSeqlevels(refGR, chrName, pruning.mode = c("coarse"))
  
  ## Output BED file
  cat("Saving BED file", "\n")
  writeGRangesToBEDwithStrand(refGR, outputPath)}

####
annotatePeakCompressor <- function(peaks){
  
  ##
  anno <- annotatePeak(
    peak = peaks,
    tssRegion = c(-1000, 500),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    level = "gene",
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
    annoDb = "org.Hs.eg.db",
    addFlankGeneInfo = FALSE,
    flankDistance = 5000,
    sameStrand = FALSE,
    ignoreOverlap = FALSE,
    ignoreUpstream = FALSE,
    ignoreDownstream = FALSE,
    overlap = "TSS",
    verbose = TRUE)
  
  ##
  return(anno)
  
}