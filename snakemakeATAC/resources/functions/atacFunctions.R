
#### IMPORTANT ####
## MUST DISABLE SCIENTIFIC NOTATION
## OR PEAKS STARTING AT COORDS ENDING IN 000 etc
## WILL BE CONVERTED TO EXPONENETS AND WILL CAUSE ERRORS
options(scipen = 999)
options(stringsAsFactors = FALSE)

######################################################################################################################################
#### Common Functions ################################################################################################################
######################################################################################################################################

#### Import a BED file as a GRanges object ####
importBED <- function(bedPath){
  suppressMessages(library(genomation))
  bedin <- readBed(bedPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
  bedin <- keepStandardChromosomes(bedin, pruning.mode = "coarse")
  bedin <- trim(bedin)
  return(bedin)}

#### Merge two or more GRanges objects ####
getGRangesUnion <- function(grList = list()){
  numItems <- length(grList)
  cat("Processing the union of", numItems, "GRanges objects", "\n")
  ## Unlist the GRanges
  for (a in 1:numItems){
    com <- paste0("tempGR", a, " <- grList[[", a, "]]")
    eval(parse(text = com))}
  ## Merge the first two
  cat("Merging first two GRanges", "\n")
  unionGR <- union(tempGR1, tempGR2)
  ## Merge the rest
  cat("Merging remaining GRanges", "\n")
  if (numItems > 2){
    for (b in 3:numItems){
      com <- paste0("unionGR <- GenomicRanges::union(unionGR, tempGR", b, ")")
      eval(parse(text = com))}}
  ## Return
  return(unionGR)}

#### Write a GRanges object to a 3-column BED file ####
writeGRangesToBED <- function(inputGR, filePath){
  
  ## BED files are 0-coordinate based, so subtract start by 1
  cat("Converting GRanges object to BED format", "\n")
  dataframeGR <- data.frame(
    chrom = seqnames(inputGR),
    chromStart = start(inputGR)-1,
    chromEnd = end(inputGR)
  )
  ## Write the output file
  cat("Writing BED output file at path:", filePath, "\n")
  write.table(
    dataframeGR,
    file = filePath,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F)}

#### Write a GRanges object to a 3-column BED file ####
writeGRangesToBEDwithStrand <- function(inputGR, filePath){
  
  ##
  numRecords <- length(inputGR)
  
  ## BED files are 0-coordinate based, so subtract start by 1
  ## Must keep name and score columns for readBed to work with strand
  cat("Converting GRanges object to BED format", "\n")
  dataframeGR <- data.frame(
    chrom = seqnames(inputGR),
    chromStart = start(inputGR)-1,
    chromEnd = end(inputGR),
    name = c(1:numRecords),
    score = rep(1, time = numRecords),
    strand = as.character(strand(inputGR)))
  ## Write the output file
  cat("Writing BED output file at path:", filePath, "\n")
  write.table(
    dataframeGR,
    file = filePath,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F)}

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

#### Retrieve underlying base sequence from a GRanges object
getUnderlyingSequence <- function(grInput){
  
  ## Report
  cat("Running getUnderlyingSequence function", "\n")
  
  ## Load required libraries
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  genome <- Hsapiens
  
  ## Get the underlying sequence
  siteSequences <- Views(genome, grInput)
  siteSequences <- as.character(siteSequences)
  
  ## Return the character sequences
  return(siteSequences)}

#### Create a fasta sequence file of hg38
generateFastaHg38 <- function(outputPath){
  
  ##
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  genome <- Hsapiens
  
  ## 
  export(genome, outputPath, format = "fasta")}

#### Convert a vector of gene symbols to entrez ids
convertSymbolsToEntrez <- function(geneSymbol){
  
  ## Report
  cat("Running convertSymbolToEntrez function", "\n")
  numSymbols <- length(geneSymbol)
  cat("Processing", numSymbols, "gene symbols", "\n")
  
  ## Load required libraries
  suppressMessages(library(org.Hs.eg.db))
  
  ## Get ENTREZ/Symbol mappings
  cat("Getting entrez/symbol mappings", "\n")
  entrezMapList <- as.list(org.Hs.egALIAS2EG)
  ## Remove entries with missing data
  entrezMapList <- entrezMapList[!is.na(entrezMapList)]
  ## Get the symbols
  allSymbols <- names(entrezMapList)
  numEntries <- length(allSymbols)
  ## Seems like the first symbol/entrez mapping is accurate, only use first id where more than one entrez id is given for a symbol
  ## Get the first id for each symbol as a vector
  firstEntrezID <- c()
  for (a in 1:numEntries){firstEntrezID[a] <- entrezMapList[[a]][1]}
  
  ## Convert the symbol
  cat("Converting symbols to entrez IDs", "\n")
  returnIDs <- c()
  for (b in 1:numSymbols){
    entrezIdx <- which(allSymbols == geneSymbol[b])
    if (length(entrezIdx) == 0){
      returnIDs[b] <- 0
    } else {
      returnIDs[b] <- firstEntrezID[entrezIdx]}}
  
  ## Return the vector
  return(returnIDs)}

#### Convert a vector of entrez IDs to gene symbols
convertEntrezToSymbol <- function(entrezID){
  
  ## Report
  cat("Running convertEntrezToSymbol function", "\n")
  numIDs <- length(entrezID)
  cat("Processing", numIDs, "entrez IDs", "\n")
  
  ## Load required libraries
  suppressMessages(library(org.Hs.eg.db))
  
  ## Get ENTREZ/Symbol mappings
  cat("Getting entrez/symbol mappings", "\n")
  entrezMapList <- as.list(org.Hs.egALIAS2EG)
  ## Remove entries with missing data
  entrezMapList <- entrezMapList[!is.na(entrezMapList)]
  ## Get the symbols
  allSymbols <- names(entrezMapList)
  numEntries <- length(allSymbols)
  ## Seems like the first symbol/entrez mapping is accurate, only use first id where more than one entrez id is given for a symbol
  ## Get the first id for each symbol as a vector
  firstEntrezID <- c()
  for (a in 1:numEntries){firstEntrezID[a] <- entrezMapList[[a]][1]}
  
  ##
  returnIDs <- c()
  
  ##
  for (b in 1:numIDs){
    
    symbolIdx <- which(firstEntrezID == entrezID[b])
    
    if (length(symbolIdx) == 0){
      
      returnIDs[b] <- 0
      
    } else {
      
      returnIDs[b] <- allSymbols[symbolIdx]}}
  
  ## Return the vector
  return(returnIDs)}

#### Retrieve a TSS window for genes *UNFINISHED*
getGeneTSSwindow <- function(entrezID = "ALL", upstream = 1000, downstream = 500){
  
  ## Report
  cat("Running getGeneLocus function", "\n")
  cat("Processing entrez IDs", entrezID, "\n")
  
  ## Load required libraries
  suppressMessages(library(Homo.sapiens))
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Get gene info and subset to standard only
  TSSinfo <- promoters(txdb, upstream = 1000, downstream = 500)
  TSSinfo <- keepStandardChromosomes(TSSinfo, pruning.mode = "coarse")
  TSSinfo <- trim(TSSinfo)
}

#### Retrieve genomic positions of genes *UNFINISHED*
getPromoterWindow <- function(entrezID = "ALL", upstream= 500, downstream = 100){
  
  ## Report
  cat("Running getGeneLocus function", "\n")
  cat("Processing entrez IDs", entrezID, "\n")
  
  ## Load required libraries
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Get gene info and subset to standard only
  geneInfo <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  geneInfo <- keepStandardChromosomes(geneInfo, pruning.mode = "coarse")
  geneInfo <- trim(geneInfo)
  geneInfo <- promoters(geneInfo, upstream = upstream, downstream = downstream)
  
  ## Return
  return(geneInfo)}

#### Split GRanges into even bins ####
splitGRangesToBin <- function(inputGR, numBins, returnBin){
  
  ## Report
  cat("Running splitGRangesToBin function", "\n")
  cat("Target number of bins", numBins, "\n")
  cat("Bin to be returned", returnBin, "\n")
  
  ## Load required libraries
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(Biostrings))
  suppressMessages(library(GenomicRanges))
  
  ## Get number of sites
  numSites <- length(inputGR)
  cat("Found", numSites, "total sites in the input GR", "\n")
  
  ## Get the indices
  grVec <- 1:numSites
  binInfo <- chunk(grVec, numBins)
  com <- paste0("binIdx <- binInfo[['", returnBin, "']]")
  eval(parse(text = com))
  binIdx <- as.numeric(binIdx)
  
  ## Add metadata indices
  score <- bindingSites@elementMetadata@listData$score
  score2 <- bindingSites@elementMetadata@listData$score2
  idx <- grVec
  df <- data.frame(score = score, score2 = score2, idx = idx)
  mcols(bindingSites, use.names = FALSE) <- df
  
  ## Subset the GR and return
  subGR <- bindingSites[(elementMetadata(bindingSites)[, "idx"] %in% binIdx)]
  return(subGR)
  
} # end splitGRangesToBin

#### Generic function for chunking an value
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

#### Get hg38 reference in GRanges
makeGRangeshg38 <- function(){
  
  ## Report
  cat("Running makeGRangeshg38 function", "\n")
  
  ## Load required libraries
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  
  ## Get the GR
  outGR <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
  
  ## Clean it up
  outGR <- keepStandardChromosomes(outGR, pruning.mode = "coarse")
  outGR <- trim(outGR)
  outGR <- sort(outGR)
  
  ## Return
  return(outGR)
}

#### Return a DF with hg38 chromosome and their lengths
getChrLengthsHg38 <- function(){
  
  chr1len <- 248956422
  chr2len <- 242193529
  chr3len <- 198295559
  chr4len <- 190214555
  chr5len <- 	181538259
  chr6len <- 170805979
  chr7len <- 159345973
  chr8len <- 145138636
  chr9len <- 	138394717
  chr10len <- 133797422
  chr11len <- 135086622
  chr12len <- 133275309
  chr13len <- 114364328
  chr14len <- 107043718
  chr15len <- 101991189
  chr16len <- 90338345
  chr17len <- 83257441
  chr18len <- 80373285
  chr19len <- 58617616
  chr20len <- 64444167
  chr21len <- 46709983
  chr22len <- 50818468
  chrXlen <- 156040895
  chrYlen <- 57227415
  #
  lengths <- c(chr1len, chr2len, chr3len, chr4len, chr5len, chr6len, chr7len, chr8len, chr9len, chr10len,
               chr11len, chr12len, chr13len, chr14len, chr15len, chr16len, chr17len, chr18len, chr19len, chr20len,
               chr21len, chr22len, chrXlen, chrYlen)
  #
  id = factor(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),
              levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                         "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  #
  df <- data.frame(chr = id, length = lengths, stringsAsFactors = FALSE)
  #
  return(df)
}

####
makeChrLabelsAsFactor <- function(){
  id = factor(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),
              levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                         "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  return(id)
}

#### Check if an ENTREZ id has a match in motifDb
checkMotifDbforEntrezID <- function(geneEntrez){
  
  ## Report
  cat("Running checkMotifDbforEntrezID function", "\n")
  cat("Input entrez ID:", geneEntrez, "\n")
  
  ## Load required libraries
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(MotifDb))
  
  ## Get ENTREZ/Symbol mappings
  cat("Getting entrez/symbol mappings", "\n")
  entrezMapList <- as.list(org.Hs.egALIAS2EG)
  ## Remove entries with missing data
  entrezMapList <- entrezMapList[!is.na(entrezMapList)]
  ## Get the symbols
  allSymbols <- names(entrezMapList)
  numEntries <- length(allSymbols)
  ## Seems like the first symbol/entrez mapping is accurate, only use first id where more than one entrez id is given for a symbol
  ## Get the first id for each symbol as a vector
  firstEntrezID <- c()
  for (a in 1:numEntries){firstEntrezID[a] <- entrezMapList[[a]][1]}
  
  ## Get all gene symbols associated with the given entrez ID
  cat("Finding all symbols mapped to current entrez ID", "\n")
  symbolIdx <- which(firstEntrezID == geneEntrez)
  numMatched <- length(symbolIdx)
  cat("Found", numMatched, "gene symbols mapped to current entrez ID", "\n")
  matchedSymbols <- allSymbols[symbolIdx]
  cat("Matching gene symbols:", matchedSymbols, "\n")
  
  ## Query the database to determine which symbols, if any, are present
  cat("Querying motifDB for matched gene symbols", "\n")
  mdbHuman <- query(MotifDb, 'hsapiens')
  foundSymbol <- c()
  ## Test each symbol
  for (b in 1:numMatched){
    geneSymbol <- matchedSymbols[b]
    geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == geneSymbol)
    if (length(geneIdx) == 0){
      foundSymbol[b] = FALSE
    } else {
      foundSymbol[b] = TRUE
    }}
  
  ## Did any of the symbols return a match?
  matchFound <- TRUE %in% foundSymbol
  return(matchFound)}

#### Clean genomic ranges, subset to standard, trim, remove mitochondrial, resort
cleanGRanges <- function(inputGRanges){
  
  ## Report
  cat("Running cleanGRanges function", "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(GenomicRanges))
  
  ## Keep only the standard chromosomes
  cleanGR <- inputGRanges
  cat("Subsetting to standard chromosomes only", "\n")
  cleanGR <- keepStandardChromosomes(cleanGR, pruning.mode = "coarse")
  
  ## Trim
  cat("Trimming", "\n")
  cleanGR <- trim(cleanGR)
  
  ## Remove chrM entries
  cat("Removing mitochondrial binding sites", "\n")
  indexChrM <- which(seqnames(cleanGR) == "chrM")
  if (length(indexChrM > 0)){cleanGR <- cleanGR[-indexChrM]}
  
  ## Resort the GRanges
  cat("Resorting seqlevels", "\n")
  cleanGR <- sortSeqlevels(cleanGR)
  cleanGR <- sort(cleanGR)
  
  ## Return the cleaned GRanges
  cat("Returning cleaned GRanges", "\n")
  return(cleanGR)}


######################################################################################################################################
#### Plotting Functions ##############################################################################################################
######################################################################################################################################

#### Generate motif insertion probability graph ####
plotFragmentSizeDistribution <- function(bamPath, outputPath1, outputPath2, outputPath3, outputPath4){
  
  ## Report
  cat("Running plotFragmentSizeDistribution function", "\n")
  cat("Bam path:", bamPath, "\n")
  cat("Output path for plot 1:", outputPath1, "\n")
  cat("Output path for plot 2:", outputPath2, "\n")
  cat("Output path for plot 3:", outputPath3, "\n")
  cat("Output path for plot 4:", outputPath4, "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(magrittr))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))
  
  ## Load the reads
  cat("Loading reads from bam file", "\n")
  atacReads <- readGAlignmentPairs(bamPath,
                                   param = ScanBamParam(mapqFilter = 1, flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE),
                                                        what = c("qname", "mapq", "isize")))
  ##
  atacReads_read1 <- GenomicAlignments::first(atacReads)
  insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
  
  ## Plot 1
  cat("Generating fragment size distribution 1", "\n")
  svg(filename = outputPath1)
  fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                              Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                                                       Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
    geom_line()
  fragLenPlot + theme_bw()
  dev.off()
  
  ## Plot 2
  cat("Generating fragment size distribution 2", "\n")
  svg(filename = outputPath2)
  fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
  dev.off()
  
  ## Plot 3 
  cat("Generating fragment size distribution 3", "\n")
  svg(filename = outputPath3)
  fragLenPlot + geom_vline(xintercept = c(180, 247),
                           colour = "red") + geom_vline(xintercept = c(315, 437),
                                                        colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
  dev.off()
  
  ## Plot 4
  cat("Generating fragment size distribution 4", "\n")
  svg(filename = outputPath4)
  fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 247),
                                                                colour = "red") + geom_vline(xintercept = c(315, 437),
                                                                                             colour = "darkblue") +
    geom_vline(xintercept = c(100),
               colour = "darkgreen") + theme_bw()
  dev.off()
}

#### Generate motif insertion probability graph ####
plotInsertionSiteProbability <- function(insertionMatrix, motifStart, motifWidth, plotTitle, svgOutPath){
  
  ## Get matrix data
  matrixTotalSignal <- sum(insertionMatrix)
  columnTotalSignal <- colSums(insertionMatrix)
  rowTotalSignal <- rowSums(insertionMatrix)
  insertionProbability <- (columnTotalSignal / matrixTotalSignal) * 100
  maxInsertionProbability <- max(insertionProbability)
  totalBP <- length(columnTotalSignal)
  flankBP <- ((totalBP - motifWidth) / 2 )
  
  ## Set the plot labels
  xLabel = "Distance to motif (bp)"
  yLabel = "Tn5 fragmentation probability"
  
  ## Set the plot axis limits
  xLimit <- c(0, totalBP + 1)
  yLimit <- c(0, maxInsertionProbability * 1.5)
  
  #### Set up the plot viewports ####
  grid.newpage()
  
  ## Set the outer bounds of the figure
  vpRegion <- plotViewport(margins=c(5.1, 5.1, 4.1, 2.1), name="plotRegion")
  pushViewport(vpRegion)
  
  ## Setup the graph plotting area
  vpGraph <- viewport(
    y = .4, 
    height = .8,
    xscale = xLimit,
    yscale = yLimit,
    name="graphRegion")
  
  pushViewport(vpGraph)
  
  ## Add the insertion probability data line
  # lwd = line width
  # col = line color
  grid.lines(
    x = 1:totalBP,
    y = insertionProbability,
    default.units = "native",
    gp = gpar(lwd = 1, col = "darkred"))
  
  ## Add the x-axis and tick marks
  # at = is a numeric vector with the x-axis locations for tick marks
  # label = specify labels on tick marks
  grid.xaxis(at = 
               c(seq(1, flankBP, length.out = 3),
                 flankBP + seq(1, motifWidth),
                 flankBP + motifWidth + seq(1, flankBP, length.out = 3)),
             label = c(-(flankBP + 1 - seq(1, flankBP + 1, length.out = 3)),
                       rep("", motifWidth),
                       seq(0, flankBP, len = 3)))
  
  ## Add the y-axis
  grid.yaxis()
  
  ## Adds the dashed lines bounding the motif
  grid.lines(
    x = c(flankBP, flankBP, 0),
    y = c(0, max(insertionProbability), yLimit[2]),
    default.units = "native",
    gp = gpar(lty = 2))
  grid.lines(
    x = c(flankBP + motifWidth + 1, flankBP + motifWidth + 1, totalBP),
    y = c(0, maxInsertionProbability, yLimit[2]),
    default.units = "native",
    gp = gpar(lty = 2))
  
  #### Add the plot title and x and y axis labels
  upViewport()
  vpTitle <- viewport(
    y = .9,
    height = .2,
    xscale = c(0, totalBP + 1),
    name = "title")
  
  pushViewport(vpTitle)
  upViewport()
  
  vpLegend <- viewport(
    x = 0.5,
    y = 0.5,
    width = convertX(unit(1, "lines"), unitTo = "npc"),
    height = convertY(unit(1, "lines"), unitTo = "npc"),
    just = c("right", "top"), name = "legendWraper")
  pushViewport(vpLegend)
  upViewport()
  
  grid.text(
    plotTitle,
    y = unit(1, "npc")-convertY(unit(1, "lines"), unitTo = "npc"),
    gp = gpar(cex = 1.2, fontface = "bold"))
  upViewport()
  
  grid.text(xLabel, y = unit(1, 'lines'))
  grid.text(yLabel, x = unit(1, 'line'), rot = 90)
  
  #### Save the image as an svg
  grid.export(svgOutPath)
  dev.off()
}

#### Generate a motif-aligned signal heatmap ####
plotMotifAlignedHeatmap <- function(insertionMatrix, bindingSites, plotTitle, svgOutPath){
  
  #### Margin controls ####
  # margin(a,b,c,d)
  # a = size of graph from top to bottom, higher value = smaller. default = 0.1
  # b = size of graph from left to right, higher value = smaller. default = 0.005
  # c = flips x axis?
  # d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
  # good settings for ATACseq = c(0.1,0.005,0.05,0.2)
  # bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values
  # upper.extreme controls to heatmap scale
  
  ## Get the matrix data
  numSites <- length(insertionMatrix[,1])
  numBP <- length(insertionMatrix[1,])
  
  ## Create list structure required to run the featureAlignedHeatmap function
  plotData <- list()
  com <- paste0("plotData$", plotTitle, " <- insertionMatrix")
  eval(parse(text = com))
  
  ## Make the plot
  ChIPpeakAnno::featureAlignedHeatmap(
    plotData,
    feature.gr = reCenterPeaks(bindingSites, width = numBP),
    upper.extreme = 1,
    annoMcols = "score",
    sortBy = "score",
    n.tile = numBP,
    margin = c(0.1, 0.005, 0.05, 0.2),
    color = colorRampPalette(c("white", "red"), bias = 0.5)(100),
    gp = gpar(fontsize = 10),
    newpage = TRUE)
  grid.export(svgOutPath)
  dev.off()
}

#### Generate peak ideogram
generatePeakIdeogram <- function(peakPath, outputPath, quantileThreshold = 0.9){
  
  ## Report
  cat("Running generatePeakIdeogram function", "\n")
  cat("Input peak file:", peakPath, "\n")
  cat("Output ideogram file:", outputPath, "\n")
  cat("Quantile threshold:", quantileThreshold, "\n")
  
  ## Load required libraries
  suppressMessages(library(RIdeogram))
  suppressMessages(library(genomation))
  suppressMessages(library(GenomicRanges))
  
  ##
  cat("Loading peaks file", "\n")
  peaks <- importBED(peakPath)
  
  ##
  cat("Determining quantiles", "\n")
  quant <- quantile(peaks@elementMetadata@listData[["score"]], probs = quantileThreshold)
  idx <- which(peaks@elementMetadata@listData[["score"]] > quant[1])
  peaks <- peaks[idx]
  
  ## Need to convert from full "chr1" values to only integer values
  cat("Converting chromosome names", "\n")
  chrNames <- c(as.character(seqnames(peaks)))
  chrNames <- gsub("chr", "", chrNames)
  
  ## Need to put data in proper format
  cat("Formatting data for plotting", "\n")
  peakData <- data.frame(
    Chr = chrNames,
    Start = peaks@ranges@start,
    End = peaks@ranges@start + peaks@ranges@width,
    Value = peaks@elementMetadata@listData[["score"]])
  
  ##
  cat("Plotting and saving ideogram", "\n")
  ideogram(karyotype = human_karyotype, overlaid = peakData, output = outputPath)
  
}

####
plotPeakAnnoCov <- function(peakPath, outputPath){
  
  ## Report
  cat("Running plotPeakAnnoCov function", "\n")
  cat("Input peak file:", peakPath, "\n")
  
  ## Load required libraries
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(ChIPseeker))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(org.Hs.eg.db))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Import the peaks
  cat("Importing peaks from BED file", "\n")
  peaks <- ChIPseeker::readPeakFile(peakPath)
  peaks <- keepStandardChromosomes(peaks)
  peaks <- trim(peaks)
  
  ## Coverage plots make plot of genome wide peak coverage
  svg(file = outputPath)
  weightname <- names(peaks@elementMetadata@listData[2])
  covplot(peaks, weightCol = weightname)
  dev.off()}

####
plotTSSpeakProfile <- function(peakPath, outputPath){
  
  ## Report
  cat("Running plotTSSpeakProfile function", "\n")
  cat("Input peak file:", peakPath, "\n")
  
  ## Load required libraries
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(ChIPseeker))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(org.Hs.eg.db))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Import the peaks
  cat("Importing peaks from BED file", "\n")
  peaks <- ChIPseeker::readPeakFile(peakPath)
  peaks <- keepStandardChromosomes(peaks)
  peaks <- trim(peaks)
  
  ## Coverage plots make plot of genome wide peak coverage
  svg(file = outputPath)
  peakHeatmap(peaks, TxDb=txdb, upstream=3000, downstream=3000, color="red")
  dev.off()}

####
plotAvgTSSpeakProfile <- function(peakPath, outputPath){
  
  ## Report
  cat("Running plotAvgTSSpeakProfile function", "\n")
  cat("Input peak file:", peakPath, "\n")
  
  ## Load required libraries
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(ChIPseeker))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(org.Hs.eg.db))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Import the peaks
  cat("Importing peaks from BED file", "\n")
  peaks <- ChIPseeker::readPeakFile(peakPath)
  peaks <- keepStandardChromosomes(peaks)
  peaks <- trim(peaks)
  
  ## Coverage plots make plot of genome wide peak coverage
  svg(file = AvgPeakProfileTSS)
  plotAvgProf2(peaks, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  dev.off()}

####
generatePeakAnnotationPlots <- function(peakPath, outputPath1, outputPath2, outputPath3, outputPath4){
  
  ## Report
  cat("Running generatePeakAnnotationPlots function", "\n")
  cat("Input peak file:", peakPath, "\n")
  
  ## Load required libraries
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(ChIPseeker))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(org.Hs.eg.db))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Import the peaks
  cat("Importing peaks from BED file", "\n")
  peaks <- ChIPseeker::readPeakFile(peakPath)
  peaks <- keepStandardChromosomes(peaks)
  peaks <- trim(peaks)
  
  ## Coverage plots make plot of genome wide peak coverage
  peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  ##
  svg(file = outputPath1)
  plotAnnoPie(peakAnno)
  dev.off()
  ##
  svg(file = outputPath2)
  plotAnnoBar(peakAnno)
  dev.off()
  ##
  svg(file = outputPath3)
  vennpie(peakAnno)
  dev.off()
  ##
  svg(file = outputPath4)
  upsetplot(peakAnno)
  dev.off()}


######################################################################################################################################
#### Binding Sites Functions #########################################################################################################
######################################################################################################################################

#### Query motifDB to return binding sites of given gene ####
getAllBindingSites <- function(geneSymbol, pwmScanScore = "95%"){
  
  ## Report
  cat("Running getAllBindingSites function", "\n")
  cat("Input gene symbol:", geneSymbol, "\n")
  cat("Input PWM matching score:", pwmScanScore, "\n")
  
  ## Load required libraries
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(Biostrings))
  suppressMessages(library(MotifDb))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(TFBSTools))
  genome <- Hsapiens
  
  #### First, try to find gene in motifDb ####
  cat("Querying motifDB", "\n")
  mdbHuman <- query(MotifDb, 'hsapiens')
  cat(length(mdbHuman), "total records found in database", "\n")
  geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == geneSymbol)
  numMatchMotifDb <- length(geneIdx)
  cat("Found", length(numMatchMotifDb), "records matching current gene", "\n")
  
  ## If found in motifDb, process
  if (numMatchMotifDb > 0){
    
    cat("Found gene in motifDB", "\n")
    
    ##
    cat("Retrieving relevant records", "\n")
    tempMotifs <- list()
    c <- 1
    for (idx in geneIdx){
      tempMotifs[c] <- mdbHuman@listData[idx]
      c <- c+1}
    
    ##
    cat("finding unique motifs", "\n")
    uniqueMotifs <- unique(tempMotifs)
    numUniqueMotifs <- length(uniqueMotifs)
    cat("found", numUniqueMotifs, "unique motifs", "\n")
    
    ##
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
    
    ## Remove chrM entries
    cat("Removing mitochondrial binding sites, if present", "\n")
    indexChrM <- which(seqnames(allSites) == "chrM")
    if (length(indexChrM > 0)){
      allSites <- allSites[-indexChrM]
    }
    
    ## Perform final trim to standard only
    cat("Trimming to standard chromosomes only", "\n")
    allSites <- keepStandardChromosomes(allSites, pruning.mode = "coarse")
    allSites <- trim(allSites)
    allSites <- sortSeqlevels(allSites)
    allSites <- sort(allSites)
    
    ##
    cat("Returning binding sites", "\n")
    return(allSites)
    
  } # end motifDB processing
  
  ## If not found in motifDb, try to load from pwm folder
  if (numMatchMotifDb == 0){
    cat("No matches found in motifDB, trying local jaspar database", "\n")
    workDir <- getwd()
    jasparPath <- paste0(workDir, "/resources/pwm/", geneSymbol, ".jaspar")
    
    if (file.exists(jasparPath)){
      cat("File found at path", jasparPath, "\n")
      cat("Importing JASPAR matrix", "\n")
      matrix <- readJASPARMatrix(jasparPath, matrixClass = "PWM")
      PWM <- matrix@listData[[1]]@profileMatrix
      ## Get binding sites and add score2
      cat("Scanning for binding sites", "\n")
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
      cat("Returning binding sites", "\n")
      return(allSites)
    } # end filecheck
  } # end if numMatchMotifDB = 0
}

#### Query motifDB to return binding sites of given gene ####
getAllBindingSitesJASPAR <- function(jasparPath, pwmScanScore = "95%"){
  
  ## Report
  cat("Running getAllBindingSitesJASPAR function", "\n")
  cat("Input jaspar filepath:", jasparPath, "\n")
  cat("Input PWM matching score:", pwmScanScore, "\n")
  
  ## Load required libraries
  suppressMessages(library(TFBSTools))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(Biostrings))
  suppressMessages(library(GenomicRanges))
  genome <- Hsapiens
  
  ## Load JASPAR file
  matrix <- readJASPARMatrix(jasparPath, matrixClass = "PWM")
  PWM <- matrix@listData[[1]]@profileMatrix
  
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
  cat("Returning binding sites", "\n")
  return(allSites)}

#### Input ENTREZ id, uses all possible gene symbols to query motifDb and return binding sites
omniGetAllBindingSites <- function(geneEntrez, pwmScanScore = "95%"){
  
  ## Report
  cat("Running omniGetAllBindingSites function", "\n")
  cat("Input entrez ID:", geneEntrez, "\n")
  cat("Input PWM matching score:", pwmScanScore, "\n")
  
  ## Load required libraries
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(Biostrings))
  suppressMessages(library(MotifDb))
  suppressMessages(library(GenomicRanges))
  genome <- Hsapiens
  
  ## Get ENTREZ/Symbol mappings
  cat("Getting entrez/symbol mappings", "\n")
  entrezMapList <- as.list(org.Hs.egALIAS2EG)
  ## Remove entries with missing data
  entrezMapList <- entrezMapList[!is.na(entrezMapList)]
  ## Get the symbols
  allSymbols <- names(entrezMapList)
  numEntries <- length(allSymbols)
  ## Seems like the first symbol/entrez mapping is accurate, only use first id where more than one entrez id is given for a symbol
  ## Get the first id for each symbol as a vector
  firstEntrezID <- c()
  for (a in 1:numEntries){firstEntrezID[a] <- entrezMapList[[a]][1]}
  
  ## Get all gene symbols associated with the given entrez ID
  cat("Finding all symbols mapped to current entrez ID", "\n")
  symbolIdx <- which(firstEntrezID == geneEntrez)
  numMatched <- length(symbolIdx)
  cat("Found", numMatched, "gene symbols mapped to current entrez ID", "\n")
  matchedSymbols <- allSymbols[symbolIdx]
  cat("Matching gene symbols:", matchedSymbols, "\n")
  
  ## Query the database to determine which symbols, if any, are present
  cat("Querying motifDB for matched gene symbols", "\n")
  mdbHuman <- query(MotifDb, 'hsapiens')
  foundSymbol <- c()
  ## Test each symbol
  for (b in 1:numMatched){
    geneSymbol <- matchedSymbols[b]
    geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == geneSymbol)
    if (length(geneIdx) == 0){
      foundSymbol[b] = FALSE
    } else {
      foundSymbol[b] = TRUE
    }}
  
  ## Did any of the symbols return a match?
  matchFound <- TRUE %in% foundSymbol
  
  ## If no matches are found, exit and return null
  if (matchFound == FALSE){
    cat("No known symbols for the current gene found in motifDb. Exiting", "\n")
    return(NULL)
  }
  
  ## If one or multiple matches are found, use the first match to get the binding sites
  if (matchFound == TRUE){
    cat("At least one symbol for current gene found in motifDb", "\n")
    matchIdx <- which(foundSymbol == TRUE)
    querySymbol <- matchedSymbols[matchIdx[1]]
    cat("Symbol of first match:", querySymbol, "\n")
    
    ##
    cat("Computing binding sites with query symbol at a PWM threshold of:", pwmScanScore, "\n")
    bindingSites <- getAllBindingSites(querySymbol, pwmScanScore)
    
    ## Return the binding sites
    cat("Returning binding sites", "\n")
    return(bindingSites)}}


######################################################################################################################################
#### Insertion Matrix Functions ######################################################################################################
######################################################################################################################################

#### Generate insertion matrix ####
generateInsertionMatrix <- function(bamFile, bindingSites, maxWidth, upstream = 100, downstream = 100){
  
  ## Turn off warnings
  options(warn=-1)
  
  ## Report
  cat("Running generateInsertionMatrix function", "\n")
  cat("Bam file:", bamFile, "\n")
  cat("Max width:", maxWidth, "\n")
  cat("Upstream search:", upstream, "\n")
  cat("Downstream search:", downstream, "\n")
  
  #### Load libraries ####
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
  genome <- Hsapiens
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
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
  return(insMatrix)}

#### Generate an insertion matrix for footprinting analysis one chromosome at a time ####
generateInsertionMatrixByChr <- function(bamFile, bindingSites, maxWidth, chrName, upstream = 100, downstream = 100){
  
  ## Turn off warnings
  options(warn=-1)
  
  ## Report
  cat("Running generateInsertionMatrix function", "\n")
  cat("Bam file:", bamFile, "\n")
  cat("Max width:", maxWidth, "\n")
  cat("Current chromosome:", chrName, "\n")
  cat("Upstream search:", upstream, "\n")
  cat("Downstream search:", downstream, "\n")
  
  #### Load libraries ####
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
  genome <- Hsapiens
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
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
  
  ## Convert to run-length encoding
  insRLE <- coverage(shiftedInsertions)
  
  ## Pull the data for the current given chromosome only
  com <- paste0("tempRLE <- list(insRLE@listData[['", chrName, "']])")
  eval(parse(text = com))
  
  insRLE@listData <- tempRLE
  extSitesList <- GRangesList(extSites)
  insViews <- Views(insRLE, extSitesList)
  
  ## Convert RLE views object to insertion matrix
  insMatrix <- as.matrix(insViews)
  
  ## Return the matrix
  return(insMatrix)}


######################################################################################################################################
#### Footprint Null Model Functions ##################################################################################################
######################################################################################################################################

#### Generate an uncorrected null footprint ####
generateNullFootprint <- function(totalSignal, bindingSite, upstream = 100, downstream = 100, iterations = 1000){
  
  ## Report
  #cat("Running generateNullFootprint function", "\n")
  #cat("Total signal:", totalSignal, "\n")
  #cat("Upstream search:", upstream, "\n")
  #cat("Downstream search:", downstream, "\n")
  #cat("Null model iterations:", iterations, "\n")
  
  ## Load required libraries
  #cat("Loading libraries", "\n")
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  
  ## Get binding site info
  #cat("Getting binding site info", "\n")
  motifWidth <- bindingSite@ranges@width
  windowLength <- motifWidth + upstream + downstream
  motifStrand <- as.character(bindingSite@strand)
  motifStart <- upstream + 1
  motifEnd <- (motifStart + motifWidth - 1)
  
  ## Create a GRanges for the extended site
  #cat("Creating GRanges object for analysis window", "\n")
  extSite <- bindingSite
  extStart <- start(extSite)
  extStartNew <- extStart - upstream
  start(extSite) <- extStartNew
  width(extSite) <- windowLength
  
  ## Generate the null model
  #cat("Generating null model", "\n")
  null <- rmultinom(n = iterations, size = totalSignal, prob = rep(1, windowLength))
  
  ## Subset for the motif region only
  motifNull <- null[motifStart:motifEnd,]
  
  ## Compute the mean signal in the motif for each null model
  #cat("Computing motif signal averages", "\n")
  avgNullMotifSignal <- c()
  for (a in 1:iterations){avgNullMotifSignal[a] <- mean(motifNull[,a])}
  
  ## Return the vector
  #cat("Returning average motif null signals", "\n")
  return(avgNullMotifSignal)
}

#### Generate a footprint null model with seqbias correction ####
generateNullFootprintSeqbiasCorrected <- function(totalSignal, bindingSite, seqbiasPath, fastaPath, upstream = 100, downstream = 100, iterations = 1000){
  
  ## Report
  #cat("Running generateNullFootprintSeqbiasCorrected function", "\n")
  #cat("Total signal:", totalSignal, "\n")
  #cat("Seqbias model path:", seqbiasPath, "\n")
  #cat("Fasta reference path:", fastaPath, "\n")
  #cat("Upstream search:", upstream, "\n")
  #cat("Downstream search:", downstream, "\n")
  #cat("Null model iterations:", iterations, "\n")
  
  ## Load required libraries
  #cat("Loading libraries", "\n")
  suppressMessages(library(seqbias))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  
  ## Get binding site info
  #cat("Getting binding site info", "\n")
  motifWidth <- bindingSite@ranges@width
  windowLength <- motifWidth + upstream + downstream
  motifStrand <- as.character(bindingSite@strand)
  motifStart <- upstream + 1
  motifEnd <- (motifStart + motifWidth - 1)
  
  ## Create a GRanges for the extended site
  #cat("Creating GRanges object for analysis window", "\n")
  extSite <- bindingSite
  extStart <- start(extSite)
  extStartNew <- extStart - upstream
  start(extSite) <- extStartNew
  width(extSite) <- windowLength
  
  ## Load the seqbias model
  #cat("Loading the seqbias model", "\n")
  seqbiasModel <- seqbias.load(fastaPath, seqbiasPath)
  
  ## Predict the insertion bias and invert to use as a probability for rmultinom
  #cat("Predicting insertion bias", "\n")
  biasPrediction <- seqbias.predict(seqbiasModel, extSite)
  biasPrediction <- biasPrediction[[1]]
  biasPrediction <- (1 / biasPrediction)
  
  ## Generate the bias corrected null models
  #cat("Generating null model", "\n")
  correctedNull <- rmultinom(n = iterations, size = totalSignal, prob = biasPrediction)
  
  ## Subset for the motif region only
  correctedMotifNull <- correctedNull[motifStart:motifEnd,]
  
  ## Compute the mean signal in the motif for each null model
  #cat("Computing motif signal averages", "\n")
  avgNullMotifSignal <- c()
  for (a in 1:iterations){avgNullMotifSignal[a] <- mean(correctedMotifNull[,a])}
  
  ## Return the vector
  #cat("Returning average motif null signals", "\n")
  return(avgNullMotifSignal)
}


######################################################################################################################################
#### Footprinting Functions ##########################################################################################################
######################################################################################################################################

#### Manually generate a single FP data object (useful if the size is not too large and can be done without splitting ins matrix)
manualGenerateFootprint <- function(geneSymbol, bamPath, peaksPath, pwmScanScore = "95%", upstream = 100, downstream = 100, outputPath){
  
  ## Report
  cat("Running manualGenerateFootprint function", "\n")
  cat("Input gene symbol:", geneSymbol, "\n")
  cat("Input bam file:", bamPath, "\n")
  cat("Input peak file:", peaksPath, "\n")
  cat("PWM matching threshold:", pwmScanScore, "\n")
  cat("Upstream search:", upstream, "\n")
  cat("Downstream search:", downstream, "\n")
  cat("Output path:", outputPath, "\n")
  
  ## Load required libraries
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
  genome <- Hsapiens
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Import peaks from file
  cat("Importing peaks", "\n")
  peaks <- importBED(peaksPath)
  
  ## Get binding sites
  cat("Finding binding sites", "\n")
  sites <- getAllBindingSites(geneSymbol, pwmScanScore)
  
  ## Subset the binding sites to peaks
  cat("Subsetting binding sites to peak regions only", "\n")
  sites <- subsetByOverlaps(sites, peaks)
  
  ##
  numSites <- length(sites@ranges)
  maxWidth <- max(sites@ranges@width)
  minWidth <- min(sites@ranges@width)
  cat("Found", numSites, "with a minimum motif width of", minWidth, "and a maximum motif width of", maxWidth, "\n")
  
  ## Generate the insertion matrix
  cat("Generating insertion matrix", "\n")
  insMatrix <- generateInsertionMatrix(bamPath, sites, maxWidth, upstream = upstream, downstream = downstream)
  
  ## Calculate the footprint statistics
  cat("Calculating footprint stats", "\n")
  fpStats <- calculateBasicFootprintStatistics(insMatrix, sites)
  
  ## Save the file
  cat("Saving file", "\n")
  save(fpStats, file = outputPath)
  
}

#### Analyze footprinting statistics given an insertion matrix ####
generateFootprintStats <- function(insMatrix, bindingSites, sampleName, geneName, upstream = 100, downstream = 100, annotationWidth = 1000){
  
  ## Prevent R from converting strings to factors (as is done in data frames)
  options(stringsAsFactors = FALSE)
  
  ## Report
  cat("Running generateFootprintStats function", "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Upstream search:", upstream, "\n")
  cat("Downstream search:", downstream, "\n")
  cat("Annotation width:", annotationWidth, "\n")
  numSites <- length(insMatrix[,1])
  cat("Total binding sites:", numSites, "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(ChIPseeker))
  suppressMessages(library(genomation))
  suppressMessages(library(stringr))
  genome <- Hsapiens
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Initiate the temporary vectors
  tempChr <- as.character(bindingSites@seqnames)
  tempStart <- as.numeric(bindingSites@ranges@start)
  tempWidth <- as.numeric(bindingSites@ranges@width)
  ##
  tempScore <- as.numeric(bindingSites@elementMetadata@listData[["score"]])
  tempScore2 <- as.numeric(bindingSites@elementMetadata@listData[["score2"]])
  ##
  tempTotalWidth <- c()
  tempTotalSignal <- c()
  tempTotalAverageSignal <- c() #
  ##
  tempTotalMotifSignal <- c() #
  tempAverageMotifSignal <- c()
  ##
  tempTotalFlankSignal <- c() #
  tempAverageFlankSignal <- c()
  ##
  tempTotalBackgroundSignal <- c() #
  tempAverageBackgroundSignal <- c()
  ##
  tempFlankAccessibility <- c() #
  tempFootprintDepth <- c() #
  tempLog2FlankAccessibility <- c() #
  tempLog2FootprintDepth <- c() #
  ##
  tempPvalue <- c()
  tempZscore <- c()
  tempBinding <- c()
  
  ## Do the binding calculations
  extSize <- as.numeric(upstream + downstream)
  flank1Start <- as.numeric(upstream - 50)
  flank1End <- as.numeric(upstream)
  background1Start <- as.numeric(1)
  background1End <- as.numeric(50)
  
  ## Loop over the binding sites
  for (b in 1:numSites){
    
    ##
    #cat("Processing site", b, "of", numSites, "\n")
    
    ##
    motifWidth <- tempWidth[b]
    motifStart <- upstream
    motifEnd <- as.numeric(motifStart + motifWidth)
    flank2Start <- as.numeric(motifEnd)
    flank2End <- as.numeric(flank2Start + 50)
    background2Start <- as.numeric(flank2End)
    background2End <- as.numeric(background2Start + 50)
    
    ##
    tempTotalWidth[b] <- as.numeric(motifWidth + extSize)
    tempTotalSignal[b] <- as.numeric(sum(insMatrix[b, 1:tempTotalWidth[b]]))
    tempTotalAverageSignal[b] <- as.numeric(tempTotalSignal[b] / tempTotalWidth[b])
    ##
    tempTotalMotifSignal[b] <- as.numeric(sum(insMatrix[b, motifStart:motifEnd])) #
    tempAverageMotifSignal[b] <- as.numeric(tempTotalMotifSignal[b] / motifWidth)
    ##
    tempTotalFlankSignal[b] <- as.numeric((sum(insMatrix[b, flank1Start:flank1End]) + sum(insMatrix[b, flank2Start:flank2End])))
    tempAverageFlankSignal[b] <- as.numeric(tempTotalFlankSignal[b] / 100)
    ##
    tempTotalBackgroundSignal[b] <- as.numeric((sum(insMatrix[b, background1Start:background1End]) + sum(insMatrix[b, background2Start:background2End])))
    tempAverageBackgroundSignal[b] <- as.numeric(tempTotalBackgroundSignal[b] / 100)
    
    ##
    if (tempTotalSignal[b] == 0){
      tempPvalue[b] <- NA
      tempZscore[b] <- NA
      tempBinding[b] <- NA}
    
    ##
    if (tempTotalSignal[b] != 0){
      
      #currentSite <- bindingSites[(elementMetadata(bindingSites)[, "idx"] %in% b)]
      #currentSite <- bindingSites[b]
      #currentSite <- bind[(elementMetadata(bind)[, "idx"] %in% 1)]
      averageNullMotifSignals <- generateNullFootprint(totalSignal = tempTotalSignal[b],
                                                       bindingSite = bindingSites[b],
                                                       upstream = upstream,
                                                       downstream = downstream,
                                                       iterations = 1000)
      averageNullMotifSignal <- mean(averageNullMotifSignals)
      motifSignals <- c(insMatrix[b, motifStart:motifEnd])
      
      ## Ttest must be in a try catch block
      ## If data are constant, will throw an error
      tryCatch({
        ttest <- t.test(motifSignals, mu = averageNullMotifSignal, alternative = "less", conf.level = 0.95)
        pvalue <- ttest$p.value
        tempPvalue[b] <- pvalue
        tempZscore[b] <- qnorm(pvalue)
        if (pvalue <= 0.05){tempBinding[b] <- 1}
        if (pvalue > 0.05){tempBinding[b] <- 0}
      }, warning = function(w) {
      }, error = function(e) {
      }, finally = {
      })
    }
  }
  
  ## Calculate the footprint scores
  tempFlankAccessibility <- as.numeric(tempAverageFlankSignal / tempAverageBackgroundSignal)
  tempLog2FlankAccessibility <- as.numeric(log2(tempFlankAccessibility))
  tempFootprintDepth <- as.numeric(tempAverageMotifSignal / tempAverageFlankSignal)
  tempLog2FootprintDepth <- as.numeric(log2(tempFootprintDepth))
  
  ## Annotate the binding sites
  cat("Annotating binding sites", "\n")
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
  
  ##
  tempAnnotation <- as.character(annotations@anno@elementMetadata@listData[["annotation"]])
  tempSymbol <- as.character(annotations@anno@elementMetadata@listData[["SYMBOL"]])
  
  ##
  sampleNameTemp <- rep(sampleName, times = numSites)
  geneNameTemp <- rep(geneName, times = numSites)
  
  ## Troubleshooting
  # cat("ID:", length(sampleNameTemp), "\n")
  # cat("chr:", length(tempChr), "\n")
  # cat("start:", length(tempStart), "\n")
  # cat("width:", length(tempWidth), "\n")
  # cat("score:", length(tempScore), "\n")
  # cat("score2:", length(tempScore2), "\n")
  # cat("totalSignal:", length(tempTotalSignal), "\n")
  # cat("averageSignal:", length(tempAverageSignal), "\n")
  # cat("averageMotifSignal:", length(tempAverageMotifSignal), "\n")
  # cat("averageFlankSignal:", length(tempAverageFlankSignal), "\n")
  # cat("averageBackgroundSignal:", length(tempAverageBackgroundSignal), "\n")
  # cat("flankAccessibility:", length(tempFlankAccessibility), "\n")
  # cat("footprintDepth:", length(tempFootprintDepth), "\n")
  # cat("binding:", length(tempBinding), "\n")
  # cat("pvalue:", length(tempPvalue), "\n")
  # cat("zscore:", length(tempZscore), "\n")
  # cat("annotation:", length(tempAnnotation), "\n")
  # cat("symbol:", length(tempSymbol), "\n")
  
  #### Convert to a dataframe before return ####
  ## This allows you to store different data types
  cat("Transferring data to dataframe", "\n")
  footprintStats <- data.frame(
    sampleID = as.character(sampleNameTemp),
    geneName = as.character(geneNameTemp),
    chr = as.character(tempChr),
    start = as.numeric(tempStart),
    width = as.numeric(tempWidth),
    score = as.numeric(tempScore),
    score2 = as.numeric(tempScore2),
    totalSignal = as.numeric(tempTotalSignal),
    averageTotalSignal = as.numeric(tempTotalAverageSignal),
    totalMotifSignal = as.numeric(tempTotalMotifSignal),
    averageMotifSignal = as.numeric(tempAverageMotifSignal),
    totalFlankSignal = as.numeric(tempTotalFlankSignal),
    averageFlankSignal = as.numeric(tempAverageFlankSignal),
    totalBackgroundSignal = as.numeric(tempTotalBackgroundSignal),
    averageBackgroundSignal = as.numeric(tempAverageBackgroundSignal),
    flankAccessibility = as.numeric(tempFlankAccessibility),
    footprintDepth = as.numeric(tempFootprintDepth),
    log2FlankAccessibility = as.numeric(tempLog2FlankAccessibility),
    log2FootprintDepth = as.numeric(tempLog2FootprintDepth),
    binding = as.numeric(tempBinding),
    pvalue = as.numeric(tempPvalue),
    zscore = as.numeric(tempZscore),
    annotation = as.character(tempAnnotation),
    symbol = as.character(tempSymbol))
  
  ##
  cat("Returning dataframe", "\n")
  return(footprintStats)}

#### Analyze footprinting statistics given an insertion matrix with seqbias correction ####
generateFootprintStatsSeqbiasCorrected <- function(insMatrix, bindingSites, sampleName, upstream = 100, downstream = 100, annotationWidth = 1000, seqbiasPath, fastaPath){
  
  ## Prevent R from converting strings to factors (as is done in data frames)
  options(stringsAsFactors = FALSE)
  
  ## Report
  cat("Running generateFootprintStats function", "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Upstream search:", upstream, "\n")
  cat("Downstream search:", downstream, "\n")
  cat("Annotation width:", annotationWidth, "\n")
  numSites <- length(insMatrix[,1])
  cat("Total binding sites:", numSites, "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(ChIPseeker))
  suppressMessages(library(genomation))
  suppressMessages(library(stringr))
  genome <- Hsapiens
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Initiate the temporary vectors
  cat("Initiating temporary vectors", "\n")
  tempChr <- as.character(bindingSites@seqnames)
  tempStart <- as.numeric(bindingSites@ranges@start)
  tempWidth <- as.numeric(bindingSites@ranges@width)
  tempScore <- as.numeric(bindingSites@elementMetadata@listData[["score"]])
  tempScore2 <- as.numeric(bindingSites@elementMetadata@listData[["score2"]])
  tempTotalWidth <- c()
  tempTotalSignal <- c()
  tempAverageSignal <- c()
  tempAverageMotifSignal <- c()
  tempAverageFlankSignal <- c()
  tempAverageBackgroundSignal <- c()
  tempPvalue <- c()
  tempZscore <- c()
  tempBinding <- c()
  
  ## Do the binding calculations
  cat("Adjusting analysis window", "\n")
  extSize <- as.numeric(upstream + downstream)
  flank1Start <- as.numeric(upstream - 50)
  flank1End <- as.numeric(upstream)
  background1Start <- as.numeric(1)
  background1End <- as.numeric(50)
  
  ## Report
  cat("extsize:", extSize, "\n")
  cat("flank1Start:", flank1Start, "\n")
  cat("flank1End:", flank1End, "\n")
  cat("background1Start:", background1Start, "\n")
  cat("background1End:", background1End, "\n")
  
  ## Loop over the binding sites
  cat("Processing binding sites", "\n")
  for (b in 1:numSites){
    
    ##
    cat("Site", b, "of", numSites, "\n")
    motifWidth <- tempWidth[b]
    motifStart <- upstream
    motifEnd <- as.numeric(upstream + motifWidth)
    flank2Start <- as.numeric(upstream + motifWidth)
    flank2End <- as.numeric(flank2Start + 50)
    background2Start <- as.numeric(flank2End)
    background2End <- as.numeric(motifWidth + extSize)
    
    ##
    tempTotalWidth[b] <- as.numeric(motifWidth + extSize)
    tempTotalSignal[b] <- as.numeric(sum(insMatrix[b, 1:tempTotalWidth[b]]))
    tempAverageSignal[b] <- as.numeric(tempTotalSignal[b] / tempTotalWidth[b])
    tempAverageMotifSignal[b] <- as.numeric(sum(insMatrix[b, (upstream:(motifEnd))] / motifWidth))
    tempAverageFlankSignal[b] <- as.numeric((sum(insMatrix[b, flank1Start:flank1End]) + sum(insMatrix[b, flank2Start:flank2End])) / 100)
    tempAverageBackgroundSignal[b] <- as.numeric((sum(insMatrix[b, background1Start:background1End]) + sum(insMatrix[b, background2Start:background2End])) / 100)
    
    ##
    if (tempTotalSignal[b] == 0){
      cat("Total signal is zero, assigning NA", "\n")
      tempPvalue[b] <- NA
      tempZscore[b] <- NA
      tempBinding[b] <- NA}
    
    if (tempTotalSignal[b] != 0){
      cat("Total signal is non-zero", "\n")
      cat("Calculating null signals", "\n")
      averageNullMotifSignals <- generateNullFootprintSeqbiasCorrected(totalSignal = tempTotalSignal[b],
                                                                       bindingSite = bindingSites[b],
                                                                       seqbiasPath = seqbiasPath,
                                                                       fastaPath = fastaPath,
                                                                       upstream = upstream,
                                                                       downstream = downstream,
                                                                       iterations = 1000)
      averageNullMotifSignal <- mean(averageNullMotifSignals)
      cat("Determining motif signal", "\n")
      motifSignals <- c(insMatrix[b, motifStart:motifEnd])
      cat("Performing t-test", "\n")
      ttest <- t.test(motifSignals, mu = averageNullMotifSignal, alternative = "less", conf.level = 0.95)
      cat("Assigning p-value and z-score", "\n")
      pvalue <- ttest$p.value
      tempPvalue[b] <- pvalue
      tempZscore[b] <- qnorm(pvalue)
      cat("Assigning binding prediction", "\n")
      if (pvalue <= 0.05){tempBinding[b] <- 1}
      if (pvalue > 0.05){tempBinding[b] <- 0}}
    
  }
  
  ## Calculate the footprint scores
  tempFlankAccessibility <- as.numeric(tempAverageFlankSignal / tempAverageBackgroundSignal)
  tempFootprintDepth <- as.numeric(tempAverageFlankSignal / tempAverageMotifSignal)
  
  ## Annotate the binding sites
  cat("Annotating binding sites", "\n")
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
  
  ##
  tempAnnotation <- as.character(annotations@anno@elementMetadata@listData[["annotation"]])
  tempSymbol <- as.character(annotations@anno@elementMetadata@listData[["SYMBOL"]])
  
  ##
  sampleNameTemp <- rep(sampleName, times = numSites)
  
  ## Troubleshooting
  cat("ID:", length(sampleNameTemp), "\n")
  cat("chr:", length(tempChr), "\n")
  cat("start:", length(tempStart), "\n")
  cat("width:", length(tempWidth), "\n")
  cat("score:", length(tempScore), "\n")
  cat("score2:", length(tempScore2), "\n")
  cat("totalSignal:", length(tempTotalSignal), "\n")
  cat("averageSignal:", length(tempAverageSignal), "\n")
  cat("averageMotifSignal:", length(tempAverageMotifSignal), "\n")
  cat("averageFlankSignal:", length(tempAverageFlankSignal), "\n")
  cat("averageBackgroundSignal:", length(tempAverageBackgroundSignal), "\n")
  cat("flankAccessibility:", length(tempFlankAccessibility), "\n")
  cat("footprintDepth:", length(tempFootprintDepth), "\n")
  cat("binding:", length(tempBinding), "\n")
  cat("pvalue:", length(tempPvalue), "\n")
  cat("zscore:", length(tempZscore), "\n")
  cat("annotation:", length(tempAnnotation), "\n")
  cat("symbol:", length(tempSymbol), "\n")
  
  #### Convert to a dataframe before return ####
  ## This allows you to store different data types
  cat("Transferring data to dataframe", "\n")
  footprintStats <- data.frame(
    sampleID = as.character(sampleNameTemp),
    chr = as.character(tempChr),
    start = as.numeric(tempStart),
    width = as.numeric(tempWidth),
    score = as.numeric(tempScore),
    score2 = as.numeric(tempScore2),
    totalSignal = as.numeric(tempTotalSignal),
    averageSignal = as.numeric(tempAverageSignal),
    averageMotifSignal = as.numeric(tempAverageMotifSignal),
    averageFlankSignal = as.numeric(tempAverageFlankSignal),
    averageBackgroundSignal = as.numeric(tempAverageBackgroundSignal),
    flankAccessibility = as.numeric(tempFlankAccessibility),
    footprintDepth = as.numeric(tempFootprintDepth),
    binding = as.numeric(tempBinding),
    pvalue = as.numeric(tempPvalue),
    zscore = as.numeric(tempZscore),
    annotation = as.character(tempAnnotation),
    symbol = as.character(tempSymbol))
  
  ##
  cat("Returning dataframe", "\n")
  return(footprintStats)
}

#### Aggregate footprint stats files into a single R list ####
aggregateFootprintStats <- function(inputFiles, geneNames){
  
  ## Prevent R from converting strings to factors (as is done in data frames)
  options(stringsAsFactors = FALSE)
  
  ## Report
  cat("Running aggregateFootprintStats function", "\n")
  numFiles <- length(inputFiles)
  cat("Number of input files:", numFiles, "\n")
  
  ##
  footprintStatsList <- list()
  
  ##
  for (a in 1:numFiles){
    
    load(inputFiles[a])
    geneName <- geneNames[a]
    ##
    com <- paste0("footprintStatsList$", gene, " <- footprintSiteStatistics")
    eval(parse(text = com))}
  
  ##
  return(footprintStatsList)
}


######################################################################################################################################
#### Seqbias Correction Functions ####################################################################################################
######################################################################################################################################

#### Generate a seqbias correction model ####
generateSeqbiasModel <- function(fastaPath, bamPath, bedPath, biasedPlotPath, correctedPlotPath, seqbiasModelPath){
  
  # refFastaPath - path to the reference fasta file
  # bamPath - path to the input bam file
  # refGR - a GRanges object containing the intervals to examine for the model
  # biasedPlotPath - output path for the biased insertion freqs
  # correctedPlotPath - output path for the corrected insertion freqs
  # biasModelPath - output path to save the correction model
  
  ## Report
  cat("Running generateSeqbiasModel function", "\n")
  cat("Reference fasta path:", fastaPath, "\n")
  cat("Bam path:", bamPath, "\n")
  cat("BED path:", bedPath, "\n")
  cat("Output path for biased plot:", biasedPlotPath, "\n")
  cat("Output path for corrected plot:", biasedPlotPath, "\n")
  cat("Output path for seqbias model:", seqbiasModelPath, "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(seqbias))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(ggplot2))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(genomation))
  
  ## Load the BED interval
  cat("Loading interval from BED file", "\n")
  refBED <- importBED(bedPath)
  # duplicate for both strands
  refBED[2] <- refBED[1]
  # set stand info
  strand(refBED) <- c("+", "-")
  refBED@strand@values <- droplevels(refBED@strand@values)
  
  ## Load the reference fasta file
  cat("Loading reference fasta", "\n")
  refFasta <- FaFile(fastaPath)
  open.FaFile(refFasta)
  refSeqIn <- scanFa(refFasta, refBED)
  
  ## Get the reverse stand indices and reverse complement
  cat("Getting negative strand indices", "\n")
  negIdx <- as.logical(refBED@strand == '-')
  refSeqIn[negIdx] <- reverseComplement(refSeqIn[negIdx])
  
  ## Count the reads in the given interval
  cat("Counting reads in given intervals", "\n")
  readCounts <- count.reads(bamPath, refBED, binary = T)
  readFreqs <- kmer.freq(refSeqIn, readCounts)
  
  ## Plot the biased sequence plot
  cat("Generating biased plot", "\n")
  svg(biasedPlotPath)
  P <- qplot(
    x = pos,
    y = freq,
    ylim = c(0, 0.5),
    color = seq,
    data = readFreqs,
    geom = "line")
  P <- P + facet_grid( seq ~ . )
  print(P)
  dev.off()
  
  ## Train the compensation model ##
  cat("Training compensation model", "\n")
  seqBiasModel <- seqbias.fit(
    fastaPath,
    bamPath,
    L = 20,
    R = 20)
  
  ## Predict sequence bias ##
  cat("Predicting insertion bias", "\n")
  ## The seqbias prediction must be given the underlying sequences you want to correct
  ## The model can be trained on chr1, etc.
  seqBiasPrediction <- seqbias.predict(
    seqBiasModel,
    refBED)
  
  ## Adjust sequence bias
  cat("Adjusting insertion counts", "\n")
  adjustedCounts <- mapply(
    FUN = `/`,
    readCounts,
    seqBiasPrediction,
    SIMPLIFY = F)
  
  ## Check and plot adjustment
  cat("Adjusting insertion frequencies", "\n")
  adjustedReadFreqs <- kmer.freq(
    refSeqIn,
    adjustedCounts)
  
  ## Generate and save the corrected freq plot
  cat("Generating corrected plot", "\n")
  svg(correctedPlotPath)
  P <- qplot(
    x = pos,
    y = freq,
    ylim = c(0.0, 0.5),
    color = seq,
    data = adjustedReadFreqs,
    geom = "line")
  P <- P + facet_grid( seq ~ . )
  print(P)
  dev.off()
  
  ## Return and save the model
  cat("Saving seqbias model", "\n")
  seqbias.save(seqBiasModel, seqbiasModelPath)}

#### Correct bias in an insertion matrix ####
correctBiasInsertionMatrix <- function(insertionMatrixData, seqbiasModelPath, refFastaPath, upstream = 100, downstream = 100){
  
  ##
  cat("Loading libraries", "\n")
  suppressMessages(library(seqbias))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  
  ## Load the seqbias model
  cat("Loading the seqbias model", "\n")
  seqbiasModel <- seqbias.load(fastaPath, seqbiasPath)
  
  ##
  cat("Loading data", "\n")
  sites <- insertionMatrixData[["bindingSites"]]
  ins <- insertionMatrixData[["insertionMatrix"]]
  numSites <- length(ins[,1])
  analysisWindow <- length(ins[1,])
  
  ## Get binding site info
  cat("Getting binding site info", "\n")
  motifStrand <- as.character(sites@strand)
  motifChr <- as.character(sites@seqnames)
  motifWidth <- max(sites@ranges@width)
  motifStart <- sites@ranges@start - upstream
  motifEnd <- motifStart + motifWidth + downstream + upstream - 1
  df <- data.frame(chr = motifChr,
                   start = motifStart,
                   end = motifEnd,
                   strand = motifStrand)
  grSite <- makeGRangesFromDataFrame(df)
  
  ##
  cat("Predicting bias", "\n")
  biasPrediction <- seqbias.predict(seqbiasModel, grSite)
  
  ## adjust the insertions
  cat("Adjusting insertions", "\n")
  adjustedInsertions <- matrix(data = NA, nrow = numSites, ncol = analysisWindow)
  for (a in 1:numSites){
    adjustedInsertions[a,] <- ins[a,] / biasPrediction[[a]]}
  
  ## return
  cat("Returning adjusted counts matrix", "\n")
  return(adjustedInsertions)
}


######################################################################################################################################
#### Gene Accessibility Functions ####################################################################################################
######################################################################################################################################

####
getGenePromoterAccessibility <- function(bamPath, entrezID = "ALL", upstream = 500, downstream = 100){
  
  ## Report
  cat("Running getGenePromoterAccessibility function", "\n")
  cat("Processing entrez IDs", entrezID, "\n")
  
  ## Load required libraries
  suppressMessages(library(Rsamtools))
  
  ## Get the promoter window for genes, resort the GR
  promoterWindows <- getPromoterWindow(upstream = upstream, downstream = downstream)
  promoterWindows <- sortSeqlevels(promoterWindows)
  promoterWindows <- sort(promoterWindows)
  numGenes <- length(promoterWindows)
  totalBP <- upstream + downstream
  
  ## Count total reads in library and calculate average signal per bp
  totalReads <- countBam(bamPath)
  totalReads <- totalReads[6]
  totalReads <- as.numeric(totalReads)
  hg38BP <- 3234830000
  avgSignalPerBP <- totalReads / hg38BP
  avgSignalPerWindow <- avgSignalPerBP * totalBP
  
  ## Scan the bam files and count reads
  params <- ScanBamParam(which = promoterWindows)
  readCounts <- countBam(bamPath, param = params)
  
  ## Add values normalized to total reads in library
  normSignal <- readCounts[,6]
  normSignal <- log2(normSignal / avgSignalPerWindow)
  readCounts <- cbind(readCounts, normSignal)
  
  ## Add the gene annotations
  geneIDs <- promoterWindows@ranges@NAMES
  geneSymbols <- convertEntrezToSymbol(geneIDs)
  readCounts <- cbind(geneIDs, readCounts)
  readCounts <- cbind(geneSymbols, readCounts)
  
  ## Return
  return(readCounts)
  
}

####
annotatePeaksToClosestGene <- function()
  
  
  ######################################################################################################################################
#### New Functions ###################################################################################################################
######################################################################################################################################

####
plotSitesDistributionChrPerMbp <- function(inputSites, geneName){
  
  ## Report
  cat("Running plotSitesDistributionChrPerMbp function", "\n")
  
  ## Load required libraries
  suppressMessages(library(ggplot2))
  suppressMessages(library(GenomicRanges))
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Get chr labels
  chrID <- makeChrLabelsAsFactor()
  
  ## Make sure the input binding sites are sorted and cleanup the GR
  names <- as.character(inputSites@seqnames)
  names(inputSites) <- names
  hg38ref <- makeGRangeshg38()
  circular <- hg38ref@seqinfo@is_circular
  gen <- hg38ref@seqinfo@genome
  inputSites@seqinfo@is_circular <- circular
  inputSites@seqinfo@genome <- gen
  #elementMetadata(inputSites) <- hg38ref@elementMetadata
  inputSites <- sort(inputSites)
  
  ## Get the number of sites per chr
  numSites <- inputSites@seqnames@lengths
  
  ## Normalize to chr length
  chrLengths <- getChrLengthsHg38() 
  numSitesNorm <- (numSites / chrLengths[,2]) * 1000000
  
  ## Convert to dataframe
  sitesDf <- data.frame(chr = chrID, numSites = numSitesNorm)
  
  ## Make title
  plotTitle <- paste0(geneName, " binding sites per Mbp")
  
  ## Make the plot
  p <- ggplot(data = sitesDf, aes(x = chr, y = numSites)) +
    geom_bar(stat = "identity", color = "blue", fill = "white", width = 1) + 
    labs(title = plotTitle)
  p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

####
getAverageFootprintStats <- function(inputStats, minSignalPerBp, upstream = 100, downstream = 100){
  
  ## Report
  cat("Running getAverageFootprintStats function", "\n")
  
  ## Load required libraries
  suppressMessages(library(GenomicRanges))
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Get the geneName
  geneName = inputStats[1,2]
  
  ## Get the extension width
  extWidth <- upstream + downstream
  
  ## Filter sites below minimum signal threshold
  filterIdx <- which(inputStats[,9] > minSignalPerBp)
  filterStats <- inputStats[filterIdx,]
  
  ## Filter out Inf and NA signals
  idxInf1 <- which(filterStats[,16] == Inf)
  idxInf2 <- which(filterStats[,17] == Inf)
  idxNA1 <- which(is.nan(filterStats[,17]))
  idxRemove <- c(idxInf1, idxInf2, idxNA1)
  filterStats <- filterStats[-idxRemove,]
  
  ## Make dataframes for bound and unbound subsets
  idxBound <- which(filterStats[,20] == 1)
  boundStats <- filterStats[idxBound,]
  unboundStats <- filterStats[-idxBound,]
  
  ## Calculate bound and unbound stats
  avgBoundSiteFlankAccessibility <- mean(boundStats[,16])
  avgBoundSiteFootprintDepth <- mean(boundStats[,17])
  avgUnboundSiteFlankAccessibility <- mean(unboundStats[,16])
  avgUnboundSiteFootprintDepth <- mean(unboundStats[,17])
  
  ##
  avgMatchScorePercentile <- mean(filterStats[,7])
  avgTotalSignal <- mean(filterStats[,8])
  avgTotalSignalPerBp <- mean(filterStats[,9])
  avgSiteMotifSignal <- mean(filterStats[,11])
  avgSiteFlankSignal <- mean(filterStats[,13])
  avgSiteBackgroundSignal <- mean(filterStats[,15])
  avgSiteFlankAccessibility <- mean(filterStats[,16])
  avgSiteFootprintDepth <- mean(filterStats[,17])
  
  ## Return the data
  returnDf <- data.frame(geneName = geneName,
                         flankAccessibility = avgSiteFlankAccessibility,
                         footprintDepth = avgSiteFootprintDepth,
                         boundFlankAccessibility = avgBoundSiteFlankAccessibility,
                         boundFootprintDepth = avgBoundSiteFootprintDepth,
                         unboundFlankAccessibility = avgUnboundSiteFlankAccessibility,
                         unboundFootprintDepth = avgUnboundSiteFootprintDepth,
                         matchScore = avgMatchScorePercentile,
                         totalSignal = avgTotalSignal,
                         avgTotalSignal = avgTotalSignalPerBp,
                         motifSignal = avgSiteMotifSignal,
                         flankSignal = avgSiteFlankSignal,
                         backgroundSignal = avgSiteBackgroundSignal)
  
  return(returnDf)
}

#### Combine average footprint stats
combineAverageFootprintStats <- function(itemNames){
  
  ## Report
  cat("Running combineAverageFootprintStats function", "\n")
  
  ## Load required libraries
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ##
  numItems <- length(itemNames)
  
  ##
  com <- paste0("combinedFPstats <- ", itemNames[1])
  eval(parse(text = com))
  
  ##
  for (a in 2:numItems){
    com <- paste0("currentFPstats <- ", itemNames[a])
    eval(parse(text = com))
    combinedFPstats <- rbind(combinedFPstats, currentFPstats)
  }
  
  ## Return
  return(combinedFPstats)
}

#### Combine average footprint stats by sample
combineAverageFootprintStatsBySamples <- function(geneName, inputNames, inputStats){
  
  ## Report
  cat("Running combineAverageFootprintStatsBySamples function", "\n")
  
  ## Load required libraries
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ##
  numSamples <- length(inputNames)
  inputNames <- data.frame(sample = inputNames)
  
  ##
  com <- paste0("combinedFPstats <- ", inputStats[1])
  eval(parse(text = com))
  
  ##
  for (a in 2:numSamples){
    com <- paste0("currentFPstats <- ", inputStats[a])
    eval(parse(text = com))
    combinedFPstats <- rbind(combinedFPstats, currentFPstats)
  }
  
  ## Add the sample names
  combinedFPstats <- cbind(inputNames, combinedFPstats)
  
  ## Return
  return(combinedFPstats)
}

####
plotAverageFootprintStats <- function(inputStats, title){
  
  ## Report
  cat("Running plotBasicFootprintStats function", "\n")
  
  ## Load required libraries
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Get gene names
  geneNames <- inputStats[,1]
  
  ##
  ggplot(inputStats, aes(x = flankAccessibility, y = footprintDepth)) +
    geom_point() + 
    xlim(0.8,1.9) +
    ylim(0.5,1.9) +
    labs(title=title) +
    geom_text_repel(label = geneNames, box.padding = 0.75)
}

####
plotAverageFootprintStatsBySamples <- function(inputStatsBySample, title, xmin = 0.8, xmax = 1.9, ymin = 0.5, ymax = 1.9){
  
  ## Report
  cat("Running plotBasicFootprintStats function", "\n")
  
  ## Load required libraries
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Get gene names
  sampleNames <- inputStatsBySample[,1]
  
  ## log2
  ggplot(inputStatsBySample, aes(x = flankAccessibility, y = footprintDepth)) +
    geom_point() + 
    xlim(xmin,xmax) +
    ylim(ymin,ymax) +
    labs(title=title) +
    geom_text_repel(label = sampleNames, box.padding = 0.75)
}

####
plotAverageFootprintStatsBySamplesBound <- function(inputStatsBySample, title, xmin = 0.8, xmax = 1.9, ymin = 0.5, ymax = 1.9){
  
  ## Report
  cat("Running plotBasicFootprintStats function", "\n")
  
  ## Load required libraries
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  
  ## Options
  options(scipen = 999)
  options(stringsAsFactors = FALSE)
  
  ## Get gene names
  sampleNames <- inputStatsBySample[,1]
  
  ## log2
  ggplot(inputStatsBySample, aes(x = boundFlankAccessibility, y = boundFootprintDepth)) +
    geom_point() + 
    xlim(xmin,xmax) +
    ylim(ymin,ymax) +
    labs(title=title) +
    geom_text_repel(label = sampleNames, box.padding = 0.75)
}

#### Generate the FP null distribution from WGS data
generateNullFPDistribution <- function(){
  
  ## Report
  cat("Running generateNullFPDistribution function", "\n")
  
  ## Load required libraries
  suppressMessages(library(GenomicRanges))
  
}

####
runCombined <- function(inputDir){
  
  inputFiles <- list.files(inputDir, full.names = TRUE)
  numFiles <- length(inputFiles)
  geneNames <- c()
  
  #
  for (a in 1:numFiles){
    load(inputFiles[a])
    geneName <- aggregatedFootprintStats[1,2]
    geneNames[a] <- geneName
    com <- paste0(geneName, "avgStats <<- getAverageFootprintStats(aggregatedFootprintStats, 0.5)")
    eval(parse(text = com))
  }
  
  ## expand item names
  itemNames <- c(paste0(geneNames, "avgStats"))
  
  ## combine
  combined <- combineAverageFootprintStats(itemNames)
  idxRemove <- which(is.nan(combined[,2]))
  combined <- combined[-idxRemove,]
  return(combined)
}

####
runCombinedGene <- function(inputDir, gene, sampleName, thresh){
  
  inputFiles <- list.files(inputDir, full.names = TRUE)
  numFiles <- length(inputFiles)
  geneNames <- c()
  
  #
  for (a in 1:numFiles){
    load(inputFiles[a])
    geneName <- aggregatedFootprintStats[1,2]
    geneNames[a] <- geneName
    com <- paste0(geneName, "avgStats <<- getAverageFootprintStats(aggregatedFootprintStats,", thresh, ")")
    eval(parse(text = com))
  }
  
  ## expand item names
  itemNames <- c(paste0(geneNames, "avgStats"))
  
  ## combine
  combined <- combineAverageFootprintStats(itemNames)
  idxRemove <- which(is.nan(combined[,2]))
  combined <- combined[-idxRemove,]
  com <- paste0("geneStats <- ", gene, "avgStats")
  eval(parse(text = com))
  geneStats[1,1] <- sampleName
  return(geneStats)
}

#### Get binding sites but with scanning function to ensure a minimum number are found
#### Uses sample peaks for subsetting
getAllBindingSitesWithMinimum <- function(geneSymbol, minimumSites, peaksPath){
  
  ## Report
  cat("Running getAllBindingSitesWithMinimum function", "\n")
  cat("Input gene symbol:", geneSymbol, "\n")
  cat("Minimum number of sites to return:", minimumSites, "\n")
  cat("Path to peaks file:", peaksPath, "\n")
  cat("Note: if minimum number of sites threshold cannot be matched, will return sites meeting 80% PWM match", "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(Biostrings))
  suppressMessages(library(MotifDb))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(genomation))
  genome <- Hsapiens
  pwmScanScore <- "80%"
  
  ## Importing accessibility peaks
  cat("Importing accessibility peaks file", "\n")
  accessibilityPeaks <- importBED(peaksPath)
  
  #### First, try to find gene in motifDb ####
  cat("Querying motifDB", "\n")
  mdbHuman <- query(MotifDb, 'hsapiens')
  geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == geneSymbol)
  numMatchMotifDb <- length(geneIdx)
  cat("Found", numMatchMotifDb, "records matching current gene", "\n")
  
  ##
  cat("Retrieving relevant records", "\n")
  tempMotifs <- list()
  c <- 1
  for (idx in geneIdx){
    tempMotifs[c] <- mdbHuman@listData[idx]
    c <- c+1}
  
  ##
  cat("finding unique motifs", "\n")
  uniqueMotifs <- unique(tempMotifs)
  numUniqueMotifs <- length(uniqueMotifs)
  cat("found", numUniqueMotifs, "unique motifs", "\n")
  
  ## If only one motif is found
  if (numUniqueMotifs == 1){
    cat("processing one unique motif", "\n")
    PWM <- uniqueMotifs[[1]]
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    allSites <- cleanGRanges(allSites)
    
    ## If more than one motif is found
  } else if (numUniqueMotifs > 1){
    cat("processing more than one unique motif", "\n")
    cat("processing motif 1", "\n")
    PWM <- uniqueMotifs[[1]]
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    for (a in 2:numUniqueMotifs){
      cat("processing motif", a, "\n")
      com <- paste0("PWM <- uniqueMotifs[[", a, "]]")
      eval(parse(text = com))
      sitesTemp <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
      allSites <- c(allSites, sitesTemp)} # end else if (numUniqueMotifs > 1)
    allSites <- cleanGRanges(allSites)} # end if (numUniqueMotifs == 1)
  
  ## Setup the while loop
  numBindingSites <- 0
  minScore <- min(allSites@elementMetadata@listData[["score"]])
  maxScore <- max(allSites@elementMetadata@listData[["score"]])
  currentScore <- maxScore
  
  while(numBindingSites <= minimumSites){
    
    ## Subset by score
    idxScore <- which(allSites@elementMetadata@listData[["score"]] >= currentScore)
    tempSites <- allSites[idxScore]
    
    ## Subset the cleaned sites 
    tempSites <- subsetByOverlaps(tempSites, accessibilityPeaks)
    numBindingSites <- length(tempSites@ranges)
    
    ##
    cat("Found:", numBindingSites, "at match score:", currentScore, "\n")
    currentScore <- (currentScore - 0.01)
    
    ## If you have reached a PWM score of 80%, break the while loop
    if (currentScore <= minScore){break}
    
  } # end while(numBindingSites < 10000)
  
  ##
  cat("Returning binding sites", "\n")
  return(tempSites)
}

#### Scan for binding sites with motifDB
scanBindingSitesMotifDB <- function(geneSymbol, PWM, pwmScanScore){
  
  ## Report
  cat("Running scanBindingSitesMotifDB function", "\n")
  cat("Input gene symbol:", geneSymbol, "\n")
  cat("Minimum PWM matching score:", pwmScanScore, "\n")
  
  ## Load required libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(Biostrings))
  suppressMessages(library(MotifDb))
  suppressMessages(library(GenomicRanges))
  genome <- Hsapiens
  
  ## Find the binding sites
  cat("Scanning for binding sites", "\n")
  bindingSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
  
  ## Add the percentile matching score
  cat("Adding score2", "\n")
  bindingSites@elementMetadata@listData$score2 <- bindingSites@elementMetadata@listData[["score"]] / max(bindingSites@elementMetadata@listData[["score"]])
  
  ## Clean the GRanges
  cat("Cleaning GRanges", "\n")
  bindingSites <- cleanGRanges(bindingSites)
  
  ## Return
  cat("Returning binding sites", "\n")
  return(bindingSites)}

