
####
findEnrichedMotifs <- function(inputPeaks, scanWidth = 200, seed = 666, maxMotifs = 2, maxWidth = 15){
  
  #### Notes ####
  {
    # Reference: https://compgenomr.github.io/book/motif-discovery.html
    # Reference: https://bioconductor.org/packages/release/bioc/html/motifRG.html
    # inputPeaks - GRanges of the regions of interest
    # scanWidth - peaks will be resized to this total size around the center of the peak
    # seed - set the seed for random shuffling of the input sequences to generate the background
    # maxMotifs - the maximum number of motifs to return
    # maxWidth - is this actually the minimum width?
    
  }
  
  #### Load libraries ####
  {
    
    cat("Loading libraries \n")
    suppressMessages(library(motifRG))
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    suppressMessages(library(TFBSTools))
    
  }
  
  #### Get input sequences
  {
    
    ##
    numSeq <- length(inputPeaks)
    cat("Detected", numSeq, "input regions \n")
    
    ##
    cat("Retrieving input sequences \n")
    set.seed(seed)
    inputSeq = getSeq(BSgenome.Hsapiens.UCSC.hg38, inputPeaks)
    
    ##
    cat("Generating background sequences \n")
    backgroundSeq = DNAStringSet(lapply(inputSeq, sample))
    
  }
  
  #### Find the enriched motifs ####
  {
    
    ##
    cat("Scanning for motifs \n")
    motifs = findMotifFgBg(fg.seq = inputSeq,
                           bg.seq = backgroundSeq,
                           max.motif = maxMotifs,
                           enriched.only = TRUE,
                           max.width = maxWidth,
                           both.strand = TRUE)
    
  }
  
  #### Refine motif edges ####
  {
    
    cat("Refining motif edges \n")
    refinedMotifs = lapply(motifs$motifs, function(x){
      motifRG::refinePWMMotifExtend(motifs = x@match$pattern, seqs = inputSeq)})
    
  }
  
  #### Return the refined motifs ####
  {
    
    cat("Returning refined motifs \n")
    return(refinedMotifs)
    
  }
  
} # end findEnrichedMotifs function

#### Make a plot of a motif ####
plotMotif <- function(inputMotif){
  
  cat("Plotting motifs \n")
  suppressMessages(library(TFBSTools))
  numMotif <- length(inputMotif)
  for(a in 1:numMotif){
    plotMotif <- inputMotif[[a]]
    seqLogo::seqLogo(plotMotif$model$prob)
  }
  
} # end plotMotif function

#### Find motif similarities ####
findMotifSimilarity <- function(inputMotif, motifID){
  
  ##
  suppressMessages(library(TFBSTools))
  suppressMessages(library(JASPAR2018))
  
  ##
  unknownMotif <- inputMotif$model$prob
  unknownPWM <- PWMatrix(ID = motifID, profileMatrix = unknownMotif)
  
  ## Retrieve motif annotation from JASPAR
  pwmLibrary = getMatrixSet(
    JASPAR2018,
    opts=list(
      collection = 'CORE',
      species    = 'Homo sapiens',
      matrixtype = 'PWM'))
  
  
  ## Compare the PWMs
  pwmSimilarity = PWMSimilarity(
    # JASPAR library
    pwmLibrary, 
    # out motif
    unknownPWM,
    # measure for comparison
    method = 'Pearson')
  
  ## Extract motif similarities
  pwm_library_list = lapply(pwmLibrary, function(x){
    data.frame(ID = ID(x), name = name(x))
  })
  
  ## combine the list into one data frame
  pwm_library_dt = dplyr::bind_rows(pwm_library_list)
  
  ## fetch the similarity of each motif to our unknown motif
  pwm_library_dt$similarity = pwmSimilarity[pwm_library_dt$ID]
  
  ## find the most similar motif in the library
  pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
  
  ## Return
  return(pwm_library_dt)
  
} # end findMotifSimilarity function
