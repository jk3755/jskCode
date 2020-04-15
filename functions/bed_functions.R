
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

#### Import a BED file as a GRanges object ####
importNarrowPeakCenterPeaks <- function(bedPath){
  
  ## Report
  cat("Running importNarrowPeakCenterPeaks function", "\n")

  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(genomation))
  
  ## 
  cat("Importing BED data from path:", bedPath, "\n")
  bedin <- readBed(bedPath, track.line = "AUTO", remove.unusual = TRUE, zero.based = TRUE)
  
  ## Center the peaks
  centeredBedin <- bedin
  newMidpoint <- floor(width(centeredBedin)/2)
  start(centeredBedin) <- start(centeredBedin) + newMidpoint
  width(centeredBedin) <- 1
  
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
  writeGRangesToBEDwithStrand(refGR, outputPath)
  
}