####
cat("Loading libraries", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(MotifDb))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(ChIPseeker))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(genomation))
suppressMessages(library(stringr))
suppressMessages(library(org.Hs.eg.db))
genome <- Hsapiens
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

####
cat("Setting vars", "\n")
geneName <- snakemake@wildcards[["gene"]]
outputDirectory <- "/ifs/scratch/c2b2/ac_lab/jk3755/atac/data/pros/lncap"
pwmScanScore = "99%"

cat("Functions", "\n")
#### Query motifDB to return binding sites of given gene ####
getAllBindingSites <- function(gene, pwmScanScore = "80%"){
  
  cat("querying motifDB", "\n")
  mdbHuman <- query(MotifDb, 'hsapiens')
  cat(length(mdbHuman), "records in database", "\n")
  geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
  cat(geneIdx, "\n")
  cat("Found", length(geneIdx), "records matching current gene", "\n")
  cat("retrieving relevant records from mdb", "\n")
  tempMotifs <- list()
  c <- 1
  for (idx in geneIdx){
    tempMotifs[c] <- mdbHuman@listData[idx]
    c <- c+1}
  cat("finding unique motifs", "\n")
  uniqueMotifs <- unique(tempMotifs)
  numUniqueMotifs <- length(uniqueMotifs)
  cat("found", numUniqueMotifs, "unique motifs", "\n")
  if (numUniqueMotifs > 1){
    cat("processing more than one unique motif", "\n")
    PWM <- uniqueMotifs[[1]]
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    allSites <- keepStandardChromosomes(allSites, pruning.mode = "coarse")
    allSites <- trim(allSites)
    allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
    for (a in 2:numUniqueMotifs){
      com <- paste0("PWM <- uniqueMotifs[[", a, "]]")
      eval(parse(text = com))
      sitesTemp <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
      sitesTemp <- keepStandardChromosomes(sitesTemp, pruning.mode = "coarse")
      sitesTemp <- trim(sitesTemp)
      sitesTemp@elementMetadata@listData$score2 <- sitesTemp@elementMetadata@listData[["score"]] / max(sitesTemp@elementMetadata@listData[["score"]])
      allSites <- c(allSites, sitesTemp)}
  } else {
    cat("processing one unique motif", "\n")
    PWM <- uniqueMotifs[[1]]
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    allSites <- keepStandardChromosomes(allSites, pruning.mode = "coarse")
    allSites <- trim(allSites)
    allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
  }
  
  ## Remove chrM entries, if present
  cat("Removing mitochondrial binding sites, if present", "\n")
  indexChrM <- which(seqnames(allSites) == "chrM")
  if (length(indexChrM > 0)){
    allSites <- allSites[-indexChrM]
  }
  
  cat("returning allSites", "\n")
  return(allSites)
}

#### Import a BED file as a GRanges object ####
importBED <- function(bedPath){
  bedin <- readBed(bedPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
  bedin <- keepStandardChromosomes(bedin, pruning.mode = "coarse")
  bedin <- trim(bedin)
  return(bedin)
}

#### Generate insertion matrix ####
generateInsertionMatrix <- function(bamFile, bindingSites, maxWidth, upstream = 100, downstream = 100){
  
  ## Expand the input binding sites to the analysis window
  extSites <- promoters(bindingSites, upstream = upstream, downstream = (downstream + maxWidth), use.names = TRUE)
  extSites <- keepStandardChromosomes(extSites, pruning.mode = "coarse")
  extSites <- trim(extSites)
  
  ## Set the parameter telling where to look in bam file
  param <- ScanBamParam(which = extSites)
  
  ## Load the relevant reads from the bam file
  bamIn <- readGAlignments(bamFile, param = param)
  
  ## Convert the GRanges and shift read ends to account for Tn5 insertion
  grIn <- granges(bamIn)
  grIn2 <- resize(grIn, width = 1)
  grPlus <- grIn2[which(strand(grIn2) == "+")]
  grMinus <- grIn2[which(strand(grIn2) == "-")]
  grPlusShifted <- shift(grPlus, shift = 4L)
  grMinusShifted <- shift(grMinus, shift = -5L)
  shiftedInsertions <- c(grPlusShifted, grMinusShifted)
  shiftedInsertions <- keepStandardChromosomes(shiftedInsertions, pruning.mode = "coarse")
  shiftedInsertions <- trim(shiftedInsertions)
  
  ## Convert to run-length encoding
  insRLE <- coverage(shiftedInsertions)
  insViews <- Views(insRLE, extSites)
  
  ## Convert RLE views object to insertion matrix
  insMatrix <- as.matrix(insViews)
  
  ## Return the matrix
  return(insMatrix)
}

#### Analyze footprinting statistics given an insertion matrix ####
calculateBasicFootprintStatistics <- function(insMatrix, bindingSites, upstream = 100, downstream = 100, annotationWidth = 1000){
  
  numSites <- length(insMatrix[,1])
  rawSiteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 19)
  
  colnames(rawSiteBasicStats) <- c("Chromosome", "motifStart", "motifWidth",
                                   "motifScore", "motifScore2",
                                   "totalSignal", "totalSignalPerBP", "motifSignalPerBP",
                                   "flankSignalPerBP", "backgroundSignalPerBP", "wideFlankSignalPerBP",
                                   "flankOverBackground", "motifOverFlank", "motifOverWideFlank",
                                   "Binding", "pValue", "zScore", "Annotation", "closestGene")
  
  rawSiteBasicStats[,1] <- as.character(bindingSites@seqnames)
  rawSiteBasicStats[,2] <- as.numeric(bindingSites@ranges@start)
  rawSiteBasicStats[,3] <- as.numeric(bindingSites@ranges@width)
  rawSiteBasicStats[,4] <- as.numeric(bindingSites@elementMetadata@listData[["score"]])
  rawSiteBasicStats[,5] <- as.numeric(bindingSites@elementMetadata@listData[["score2"]])
  
  ## Annotate
  annotations <- annotatePeak(
    peak = bindingSites,
    tssRegion = c(-annotationWidth, annotationWidth),
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
  
  rawSiteBasicStats[,18] <- annotations@anno@elementMetadata@listData[["annotation"]]
  rawSiteBasicStats[,19] <- annotations@anno@elementMetadata@listData[["SYMBOL"]]
  
  ####
  for (b in 1:numSites){
    
    ##
    cat("Processing site", b, "of", numSites, "\n")
    
    tempWidth <- as.numeric(rawSiteBasicStats[b,3])
    tempTotalWidth <- as.numeric(tempWidth + upstream + downstream)
    tempTotalSignal <- as.numeric(sum(insMatrix[b,1:tempTotalWidth]))
    tempTotalSignalPerBP <- as.numeric(tempTotalSignal/tempTotalWidth)
    tempMotifSignalPerBP <- as.numeric(sum(insMatrix[b,(upstream:(upstream + tempWidth))] / tempWidth))
    tempFlankSignalPerBP <- as.numeric((sum(insMatrix[b,(upstream-50):upstream])+sum(insMatrix[b,(upstream+tempWidth):((upstream+50)+tempWidth)])) / 100)
    tempBackgroundSignalPerBP <- as.numeric((sum(insMatrix[b,1:50])+sum(insMatrix[b,((upstream+downstream)+tempWidth-50):((upstream+downstream)+tempWidth)])) / 100)
    tempWideFlankSignalPerBP <- as.numeric((sum(insMatrix[b,1:upstream])+sum(insMatrix[b,(upstream+tempWidth):((upstream+downstream)+tempWidth)])) / (upstream+downstream))
    ##
    rawSiteBasicStats[b,6] <- as.numeric(tempTotalSignal)
    rawSiteBasicStats[b,7] <- as.numeric(tempTotalSignalPerBP)
    rawSiteBasicStats[b,8] <- as.numeric(tempMotifSignalPerBP)
    rawSiteBasicStats[b,9] <- as.numeric(tempFlankSignalPerBP)
    rawSiteBasicStats[b,10] <- as.numeric(tempBackgroundSignalPerBP)
    rawSiteBasicStats[b,11] <- as.numeric(tempWideFlankSignalPerBP)
    rawSiteBasicStats[b,12] <- as.numeric(tempFlankSignalPerBP / tempBackgroundSignalPerBP)
    rawSiteBasicStats[b,13] <- as.numeric(tempMotifSignalPerBP / tempFlankSignalPerBP)
    rawSiteBasicStats[b,14] <- as.numeric(tempMotifSignalPerBP / tempWideFlankSignalPerBP)
    
    if (tempTotalSignal == 0){
      
    } else {
      averageNullMotifSignal <- generateNullFP(1000, tempTotalSignal, tempTotalWidth, tempWidth)
      motifSignals <- c(insMatrix[b,(upstream:(upstream + tempWidth))])
      
      result = tryCatch({
        ttest <- t.test(motifSignals, mu = averageNullMotifSignal, alternative = "less", conf.level = 0.95)
        pvalue <- ttest$p.value
        rawSiteBasicStats[b,16] <- pvalue
        rawSiteBasicStats[b,17] <- qnorm(pvalue)
        if (pvalue < 0.05){rawSiteBasicStats[b,15] <- 1} else {rawSiteBasicStats[b,15] <- 0}
      }, warning = function(w) {
      }, error = function(e) {
      }, finally = {
      })}
  }
  
  #### Convert to a dataframe before return ####
  ## This allows you to store different data types
  cat("Transferring data to dataframe", "\n")
  rawFootprintStats <- data.frame(
    chr = as.character(rawSiteBasicStats[,1]),
    start = as.numeric(rawSiteBasicStats[,2]),
    width = as.numeric(rawSiteBasicStats[,3]),
    score = as.numeric(rawSiteBasicStats[,4]),
    score2 = as.numeric(rawSiteBasicStats[,5]),
    totalSignal = as.numeric(rawSiteBasicStats[,6]),
    totalSignalPerBP = as.numeric(rawSiteBasicStats[,7]),
    motifSignalPerBP = as.numeric(rawSiteBasicStats[,8]),
    flankSignalPerBP = as.numeric(rawSiteBasicStats[,9]),
    backgroundSignalPerBP = as.numeric(rawSiteBasicStats[,10]),
    wideFlankSignalPerBP = as.numeric(rawSiteBasicStats[,11]),
    flankAccessibility = as.numeric(rawSiteBasicStats[,12]),
    footprintDepth = as.numeric(rawSiteBasicStats[,13]),
    wideAccessibility = as.numeric(rawSiteBasicStats[,14]),
    binding = as.numeric(rawSiteBasicStats[,15]),
    pvalue = as.numeric(rawSiteBasicStats[,16]),
    zscore = as.numeric(rawSiteBasicStats[,17]),
    annotation = as.character(rawSiteBasicStats[,18]),
    symbol = as.character(rawSiteBasicStats[,19])
  )
  cat("Returning dataframe", "\n")
  return(rawFootprintStats)
}

#### Generate a null model of a TF binding footprint for hypothesis testing ####
generateNullFP <- function(iterations, inputSignal, analysisWidth, motifWidth){
  averages <- c()
  for (a in 1:iterations){
    null <- c(1:(analysisWidth))
    null <- c(as.vector(rmultinom(1, size=inputSignal, prob=rep(1, length(null)))))
    motifStart <- ((analysisWidth - motifWidth)/2)
    motifEnd <- (motifStart + motifWidth)
    motifAvg <- (sum(null[motifStart:motifEnd])) / motifWidth
    averages[a] <- motifAvg
  }
  return(mean(averages))
}

####
runGeneListFP <- function(symbols, pwmScanScore = "99%", outputDirectory){
  
  #### Bam paths ####
  LNCaP.CR01.bam.path <- snakemake@input[[2]]
  LNCaP.CR02.bam.path <- snakemake@input[[3]]
  LNCaP.CR04.bam.path <- snakemake@input[[4]]
  LNCaP.CR05.bam.path <- snakemake@input[[5]]
  LNCaP.CR07.bam.path <- snakemake@input[[6]]
  LNCaP.CR08.bam.path <- snakemake@input[[7]]
  LNCaP.WT01.bam.path <- snakemake@input[[8]]
  LNCaP.WT02.bam.path <- snakemake@input[[9]]
  
  #### Peaks ####
  peaksPath <- snakemake@input[[1]]
  LNCaP.peaks <- importBED(peaksPath)
  
  ####
  outputPath <- paste0(outputDirectory, "/fpdata.", geneName, ".RDS")
  cat("Processing", geneName, "\n")
  cat("Output path", outputPath, "\n")
  
  tryCatch({

    ## Get the binding sites and subset by peaks
    sites <- getAllBindingSites(geneName, pwmScanScore)
    sites <- subsetByOverlaps(sites, LNCaP.peaks)
    numSites <- length(sites@ranges)
    maxWidth <- max(sites@ranges@width)
    minWidth <- min(sites@ranges@width)
    
    ## Footprint stats
    LNCaP.CR01.insMatrix <- generateInsertionMatrix(LNCaP.CR01.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.CR01.footprintStats <- calculateBasicFootprintStatistics(LNCaP.CR01.insMatrix, sites)
    LNCaP.CR02.insMatrix <- generateInsertionMatrix(LNCaP.CR02.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.CR02.footprintStats <- calculateBasicFootprintStatistics(LNCaP.CR02.insMatrix, sites)
    LNCaP.CR04.insMatrix <- generateInsertionMatrix(LNCaP.CR04.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.CR04.footprintStats <- calculateBasicFootprintStatistics(LNCaP.CR04.insMatrix, sites)
    LNCaP.CR05.insMatrix <- generateInsertionMatrix(LNCaP.CR05.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.CR05.footprintStats <- calculateBasicFootprintStatistics(LNCaP.CR05.insMatrix, sites)
    LNCaP.CR07.insMatrix <- generateInsertionMatrix(LNCaP.CR07.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.CR07.footprintStats <- calculateBasicFootprintStatistics(LNCaP.CR07.insMatrix, sites)
    LNCaP.CR08.insMatrix <- generateInsertionMatrix(LNCaP.CR08.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.CR08.footprintStats <- calculateBasicFootprintStatistics(LNCaP.CR08.insMatrix, sites)
    LNCaP.WT01.insMatrix <- generateInsertionMatrix(LNCaP.WT01.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.WT01.footprintStats <- calculateBasicFootprintStatistics(LNCaP.WT01.insMatrix, sites)
    LNCaP.WT02.insMatrix <- generateInsertionMatrix(LNCaP.WT02.bam.path, sites, maxWidth, upstream = 100, downstream = 100)
    LNCaP.WT02.footprintStats <- calculateBasicFootprintStatistics(LNCaP.WT02.insMatrix, sites)
    
    ## Put in a list
    
    LNCaP.FPstats <- list()
    LNCaP.FPstats$WT01 <- LNCaP.WT01.footprintStats
    LNCaP.FPstats$WT02 <- LNCaP.WT02.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats
    LNCaP.FPstats$CR01 <- LNCaP.CR01.footprintStats

    ## save
    saveRDS(LNCaP.FPstats, file = outputPath)
  }, finally = {
  }) 
}

## run
cat("Running", "\n")
runGeneListFP(geneName, outputDirectory = outputDirectory)