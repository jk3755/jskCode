
#### Import bigwig file from disk ####
importBigwig <- function(filepath, clean = TRUE){
  
  ## Load libraries
  suppressMessages(library(rtracklayer))
  
  ## Import
  temp <- rtracklayer::import.bw(filepath)
  
  ## Clean if true
  if (clean == TRUE){temp <- cleanGRanges(temp)}
  
  ## Return
  return(temp)
  
} # end importBigwig function

#### Write bigwig to file ####
writeBigwigToFile <- function(bigwig, filepath){
  
  ## Load libraries
  suppressMessages(library(rtracklayer))
  
  ##
  rtracklayer::export.bw(bigwig, filepath)
  
} # end exportBigwig function

#### Generate averaged bigwig from multiple input bigwigs and save to disk ####
averageMultipleBigwigs <- function(bigwigList, binSize, validate = TRUE, verbose = TRUE){
  
  #### Reporting ####
  {
    cat("NOTE: it is CRITICAL that the bin size is set to the same bin size used when generating the input bigwig files \n")
  }
  
  #### Setup ####
  {
    
    ##
    if(verbose){cat("Loading libraries", "\n")}
    suppressMessages(library(GenomicRanges))
    
    ##
    numItems <- length(bigwigList)
    if(verbose){cat("Detected", numItems, "input samples", "\n")}
    
    ##
    if(verbose){cat("Bin size is set at", binSize, "bp", "\n")}
    
    ##
    if(verbose){cat("Generating genomic bins \n")}
    binRanges <- generateGenomicWindows(binSize)
    numBins <- length(binRanges)
    if(verbose){cat("Genome split to", numBins, "bins",  "\n")}
    
  }
  
  #### Initialize the data matrix ####
  {
    
    ## Initiate a matrix to store the data
    ## One row for each sample, one column for each bin
    if(verbose){cat("Initializing data matrix \n")}
    binData <- matrix(data = NA, nrow = numItems, ncol = numBins)
    
  }
  
  #### Loop over each sample and copy the bin scores to the data matrix ####
  {
    
    ##
    if(verbose){cat("Copying sample bin data to data matrix \n")}
    
    ##
    for (a in 1:numItems){
      
      ##
      if(verbose){cat("Processing bigwig", a, "\n")}
      currentBigwig <- bigwigList[[a]]
      currentIdxs <- findOverlaps(binRanges, currentBigwig)
      currentIdxs <- currentIdxs@to
      binData[a,] <- currentBigwig@elementMetadata@listData[["score"]][currentIdxs]
      
    } # end for (a in 1:numItems)
    
  }
  
  #### Find the average value of each bin ####
  {
    
    ##
    if(verbose){cat("Determining average score in each bin \n")}
    avgScores <- colMeans(binData)
    
  }
  
  #### Return the averaged bigwig object ####
  {
    
    ##
    if(verbose){cat("Generating averaged bigwig object and returning \n")}
    mcols(binRanges)$score <- avgScores
    return(binRanges)
    
  }
  
} # end averageMultipleBigwigs function
