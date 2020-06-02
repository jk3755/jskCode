

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
  return(returnIDs)
  } # end convertSymbolsToEntrez

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
  cat("Running getGeneTSSwindow function", "\n")
  cat("Processing entrez IDs", entrezID, "\n")
  
  ## Load required libraries
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Get gene info and subset to standard only
  TSSinfo <- promoters(txdb, upstream = upstream, downstream = downstream)
  
  ## Clean
  TSSinfo <- cleanGRangesBasic(TSSinfo)
  
  ## Return
  return(TSSinfo)
  
} # end getGeneTSSwindow

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
  return(geneInfo)
} # end getPromoterWindow