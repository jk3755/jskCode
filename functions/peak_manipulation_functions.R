
####
mergeAndAverageReplicatedPeaks <- function(peaks1, peaks2, minOverlap, extendPeak, verbose = FALSE){
  
  ## peak1, peak2: input peaks as GRanges, must have score value
  ## minOverlap: minimum bp overlap to be considered a hit
  ## extendPeak: bp to extend the centered merged peak up and downstream
  
  ##
  overlaps <- findOverlaps(peaks1,
                           peaks2, 
                           minoverlap = minOverlap)
  
  ##
  numOverlaps <- length(overlaps)
  cat("Detected", numOverlaps, "overlapping peaks at current threshold \n")
  
  ## Initialize the vectors to hold merged peaks data
  vecMergedChr <- c()
  vecMergedStart <- c()
  vecMergedWidth <- c()
  vecMergedScore <- c()
  
  ## Overlap indices
  idxOverlap1 <- overlaps@from
  idxOverlap2 <- overlaps@to
  
  ## For each overlap, pull both corresponding ranges
  ## Find the middle distance between their centers
  ## Set middle as new peak midpoint
  ## Extend peak upstream and downstream
  for(a in 1:numOverlaps){
    
    ##
    if(verbose){cat("Processing peak set", a, "of", numOverlaps, "total overlaps \n")}
    
    ## Idx
    idx1 <- idxOverlap1[a]
    idx2 <- idxOverlap2[a]
    
    ##
    currentPeak1 <- peaks1[idx1]
    currentPeak2 <- peaks2[idx2]
    
    ## verify that the chromosomes match
    chr1 <- as.character(seqnames(currentPeak1))
    chr2 <- as.character(seqnames(currentPeak2))
    stopifnot(chr1==chr2)
    mergedChr <- chr1
    
    ##
    start1 <- currentPeak1@ranges@start
    start2 <- currentPeak2@ranges@start
    
    ##
    peakStartMidpoint <- floor((start1 + start2) / 2)
    
    ##
    mergedStart <- peakStartMidpoint - (floor(extendPeak/2))
    mergedWidth <- (extendPeak - 1)
    
    ## Scores
    score1 <- currentPeak1@elementMetadata@listData[["score"]]
    score2 <- currentPeak2@elementMetadata@listData[["score"]]
    
    ## merged score
    mergedScore <- floor((score1 + score2) / 2)
    
    ## Update the vectors
    vecMergedChr[a] <- mergedChr
    vecMergedStart[a] <- mergedStart
    vecMergedWidth[a] <- mergedWidth
    vecMergedScore[a] <- mergedScore
    
  }
  
  #### Make the dataframe
  dfMergedPeaks <- data.frame(chromosome = vecMergedChr,
                              start = vecMergedStart,
                              end = (vecMergedStart + vecMergedWidth),
                              score = vecMergedScore)
  
  #### Convert to GRanges
  mergedPeaks <- makeGRangesFromDataFrame(dfMergedPeaks, keep.extra.columns = TRUE)
  
  #### Return
  return(mergedPeaks)

} # end mergeAndAverageReplicatedPeaks function