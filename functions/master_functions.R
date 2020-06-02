
#############################################################################################################################################################################
#### ARACNe and VIPER #######################################################################################################################################################
#############################################################################################################################################################################
{
  
  ####
  convert_regulon_from_entrez_to_symbol <- function(regulon){
    
    ## Load libraries
    cat("Loading libraries \n")
    suppressMessages(library(org.Hs.eg.db))
    
    ## Get mappings
    #cat("Getting entrez/symbol mappings \n")
    
    entrezMapList <- as.list(org.Hs.egALIAS2EG)
    ## Remove entries with missing data
    entrezMapList <- na.omit(entrezMapList)
    #entrezMapList <- entrezMapList[!is.na(entrezMapList)]
    ## Get the symbols
    allSymbols <- names(entrezMapList)
    numEntries <- length(allSymbols)
    ## Seems like the first symbol/entrez mapping is accurate, only use first id where more than one entrez id is given for a symbol
    ## Get the first id for each symbol as a vector
    firstEntrezID <- c()
    for (a in 1:numEntries){firstEntrezID[a] <- entrezMapList[[a]][1]}
    ## Put the mappings into a dataframe
    mapDF <- data.frame(entrez = firstEntrezID,
                        symbol = allSymbols,
                        stringsAsFactors = FALSE)
    
    ## Convert the regulator IDs
    #cat("Converting regulator IDs to symbol \n")
    regulatorEntrez <- names(regulon)
    regulatorSymbol <- c()
    
    ##
    for(a in 1:length(regulatorEntrez))
    {
      symbolIdx <- which(mapDF[,1] == regulatorEntrez[a])
      if (length(symbolIdx) == 0)
      {regulatorSymbol[a] <- "NULL"
      } else {regulatorSymbol[a] <- mapDF[symbolIdx[1],2]}
    }
    
    ##
    names(regulon) <- regulatorSymbol
    
    ## Convert the regulator targets
    #cat("Converting regulator target IDs to symbols \n")
    numRegulators <- length(regulon)
    ##
    for(b in 1:numRegulators){
      ##
      #cat("Processing regulator", b, "\n")
      ##
      tempRegulon <- regulon[b]
      tempTargets <- names(tempRegulon[[1]][["tfmode"]])
      targetSymbol <- c()
      ##
      for(c in 1:length(tempTargets))
      {
        symbolIdx <- which(mapDF[,1] == tempTargets[c])
        if (length(symbolIdx) == 0)
        {targetSymbol[c] <- "NULL"
        } else {targetSymbol[c] <- mapDF[symbolIdx[1],2]}
      }
      ##
      names(tempRegulon[[1]][["tfmode"]]) <- targetSymbol
      regulon[b] <- tempRegulon
    }
    #
    ## Return the converted regulon
    cat("Returning converted regulon \n")
    return(regulon)
    
  } # end convert_regulon_from_entrez_to_symbol
  
  ####
  convert_regulon_from_ensembl_to_symbol <- function(regulon, verbose = FALSE){
    
    ## load libraries
    suppressMessages(library(org.Hs.eg.db))
    
    ## get regulator names
    regulator_names <- names(regulon)
    
    ## map to symbols
    regulator_symbols <- convert_ensembl_to_symbol(regulator_names)
    
    ## convert the regulator names
    names(regulon) <- regulator_symbols
    
    ## Convert the regulator targets
    num_regulators <- length(regulon)
    
    ##
    for(a in 1:num_regulators){
      
      ##
      if(verbose){cat("Processing regulator", a, "of", num_regulators, "\n")}
      
      ##
      temp_regulon <- regulon[a]
      temp_targets <- names(temp_regulon[[1]][["tfmode"]])
      
      ##
      target_symbols <- convert_ensembl_to_symbol(temp_targets)
      
      ##
      names(temp_regulon[[1]][["tfmode"]]) <- target_symbols
      
      ##
      regulon[a] <- temp_regulon}
    
    ##
    return(regulon)
    
  } # end convert_regulon_from_entrez_to_symbol
  
  ####
  convert_regulon_from_ensembl_to_entrez <- function(regulon, verbose = FALSE){
    
    ## load libraries
    suppressMessages(library(org.Hs.eg.db))
    
    ## get regulator names
    regulator_names <- names(regulon)
    
    ## map to symbols
    regulator_symbols <- convert_ensembl_to_entrez(regulator_names)
    
    ## convert the regulator names
    names(regulon) <- regulator_symbols
    
    ## Convert the regulator targets
    num_regulators <- length(regulon)
    
    ##
    for(a in 1:num_regulators){
      
      ##
      if(verbose){cat("Processing regulator", a, "of", num_regulators, "\n")}
      
      ##
      temp_regulon <- regulon[a]
      temp_targets <- names(temp_regulon[[1]][["tfmode"]])
      
      ##
      target_symbols <- convert_ensembl_to_entrez(temp_targets)
      
      ##
      names(temp_regulon[[1]][["tfmode"]]) <- target_symbols
      
      ##
      regulon[a] <- temp_regulon}
    
    ##
    return(regulon)
    
  } # end convert_regulon_from_ensembl_to_entrez
  
}
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#############################################################################################################################################################################
#### ID conversion ##########################################################################################################################################################
#############################################################################################################################################################################
{
  
  ####
  convert_symbol_to_entrez <- function(symbols){
    
    ##
    suppressMessages(library(org.Hs.eg.db))
    
    ##
    db <- org.Hs.eg.db
    
    ##
    mapping <- mapIds( db , keys = symbols , column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    
    ##
    x <- mapping[ match( symbols , names(mapping) ) ]
    x <- as.character(x)
    x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
    
    ##
    return(x)
    
  } # end convert_symbol_to_entrez
  
  ####
  convert_entrez_to_symbol <- function(entrez){
    
    ##
    suppressMessages(library(org.Hs.eg.db))
    
    ##
    db <- org.Hs.eg.db
    
    ##
    mapping <- mapIds( db , keys = entrez , column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    
    ##
    x <- mapping[ match( entrez , names(mapping) ) ]
    x <- as.character(x)
    x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
    
    ##
    return(x)
    
  } # end convert_entrez_to_symbol
  
  ####
  convert_ensembl_to_symbol <- function(ensembl){
    
    ##
    suppressMessages(library(org.Hs.eg.db))
    
    ##
    db <- org.Hs.eg.db
    
    ##
    mapping <- mapIds( db , keys = ensembl , column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    
    ##
    #message("- Found NA - " , sum(is.na(mapping)) , " total" )
    x <- mapping[ match( ensembl , names(mapping) ) ]
    x <- as.character(x)
    x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
    
    ##
    return(x)
    
  } # end convert_ensembl_to_symbol
  
  ####
  convert_ensembl_to_entrez <- function(ensembl){
    
    ##
    suppressMessages(library(org.Hs.eg.db))
    
    ##
    db <- org.Hs.eg.db
    
    ##
    mapping <- mapIds( db , keys = ensembl , column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
    
    ##
    message("- Found NA - " , sum(is.na(mapping)) , " total" )
    x <- mapping[ match( ensembl , names(mapping) ) ]
    x <- as.character(x)
    x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
    
    ##
    return(x)
    
  } # end convert_ensembl_to_entrez
  
}
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#############################################################################################################################################################################
#### BED files ##############################################################################################################################################################
#############################################################################################################################################################################
{
  
  ####
  import_narrowpeak_as_granges_standardized <- function(filepath, width = 500, remove_overlaps = TRUE){
    
    ## Load libraries
    cat("Loading libraries \n")
    suppressMessages(library(chromVAR))
    suppressMessages(library(GenomicRanges))
    
    ## Import the peaks using chromVar function
    cat("Reading narrowPeak files \n")
    peaks <- readNarrowpeaks(filepath, width = width, non_overlapping = remove_overlaps)
    
    ## Clean
    cat("Cleaning GenomicRanges peak object \n")
    peaks <- clean_granges(peaks)
    
    ## Return
    return(peaks)
    
  } # end import_narrowpeak_as_granges_standardized
  
}
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#############################################################################################################################################################################
#### GenomicRanges ##########################################################################################################################################################
#############################################################################################################################################################################
{
  
  #### Subset a GRanges to standard chromosomes, remove chrM, trim, resort ####
  clean_granges <- function(ranges, rm.chrM = TRUE, rm.chrY = TRUE, rm.chrX = FALSE){
    
    ##
    suppressMessages(library(GenomicRanges))
    
    ## Subset to standard chromosomes only
    ranges <- keepStandardChromosomes(ranges, pruning.mode = "coarse")
    
    ## Drop the ChrM seqlev
    if (rm.chrM == TRUE){ranges <- dropSeqlevels(ranges, "chrM", pruning.mode = "coarse")}
    
    ## Drop the ChrY seqlev
    if (rm.chrY == TRUE){ranges <- dropSeqlevels(ranges, "chrY", pruning.mode = "coarse")}
    
    ## Drop the ChrX seqlev
    if (rm.chrX == TRUE){ranges <- dropSeqlevels(ranges, "chrX", pruning.mode = "coarse")}
    
    ## Trim
    ranges <- trim(ranges)
    
    ## Resort
    ranges <- sortSeqlevels(ranges)
    ranges <- sort(ranges)
    
    ##
    return(ranges)
    
  } # end clean_granges
  
  ####
  find_and_average_replicated_granges <- function(ranges1, ranges2, min_overlap, extend, verbose = FALSE){
    
    ##
    overlaps <- findOverlaps(ranges1,
                             ranges2, 
                             minoverlap = min_overlap)
    
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
      currentPeak1 <- ranges1[idx1]
      currentPeak2 <- ranges2[idx2]
      
      ## verify that the chromosomes match
      chr1 <- as.character(seqnames(currentPeak1))
      chr2 <- as.character(seqnames(currentPeak2))
      stopifnot(chr1==chr2)
      mergedChr <- chr1
      
      ##
      start1 <- currentPeak1@ranges@start
      start2 <- currentPeak2@ranges@start
      
      ##
      mergedStart <- floor((start1 + start2) / 2)
      mergedWidth <- (extend - 1)
      
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
    
  } # end find_and_average_replicated_granges
  
  ####
  get_granges_union <- function(rangesList = list()){
    
    ##
    suppressMessages(library(GenomicRanges))
    
    ##
    temp.merged <- rangesList[[1]]
    numItems <- length(rangesList)
    
    ##
    for (a in 2:numItems){
      temp.gr <- rangesList[[a]]
      temp.merged <- GenomicRanges::union(temp.merged, temp.gr)
    }
    
    ##
    return(temp.merged)
    
  } # end get_granges_union
  
  ####
  split_granges_to_promoter_and_distal <- function(ranges, promoter_width = 1000){
    
    ## a distance of 1000 bp upstream is the approx max size
    ## though individual promoters can vary from 100-1000 bp
    
    ##
    suppressMessages(library(GenomicRanges))
    
    ##
    promoterWindows <- get_gene_promoter_windows(size = promoter_width)
    
    ##
    overlaps <- mergeByOverlaps(ranges, promoterWindows, ignore.strand = TRUE)
    promoter_ranges <- overlaps@listData[["ranges"]]
    gene_symbols <- overlaps@listData[["promoterWindows"]]@elementMetadata@listData[["hgnc_symbol"]]
    entrez_ids <- overlaps@listData[["promoterWindows"]]@elementMetadata@listData[["entrez_id"]]
    mcols(promoter_ranges)$symbol <- gene_symbols
    mcols(promoter_ranges)$entrez <- entrez_ids
    promoter_ranges <- unique(promoter_ranges)
    distal_ranges <- setdiff(ranges, promoter_ranges)
    length(promoter_ranges) + length(distal_ranges)
    length(ranges)
    
    ## Convert to list
    peaksList <- list()
    peaksList$promoterPeaks <- promoter_ranges
    peaksList$distalPeaks <- distal_ranges
    
    ## Return
    return(peaksList)
    
  } # end split_granges_to_promoter_and_distal
  
  ####
  get_gene_promoter_windows <- function(size = 1000){
    
    gene_ranges <- generate_all_gene_annotation_as_granges()
    gene_promoters <- promoters(gene_ranges, upstream = size, downstream = 0)
    return(gene_ranges)
    
  } # end get_gene_promoter_windows
  
  ####
  generate_all_gene_annotation_as_granges <- function(){
    
    ##
    suppressMessages(library(biomaRt))
    suppressMessages(library(GenomicRanges))
    options(stringsAsFactors = FALSE)
    
    ##
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    attributes <- c("entrezgene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand")
    
    ##
    info <- getBM(attributes = attributes, mart = mart)
    chrNames <- c(1:22, "X", "Y")
    idx <- which(info$chromosome_name %in% chrNames)
    info_standard <- info[idx,]
    idx <- which(info_standard$hgnc_symbol == "")
    info_standard <- info_standard[-idx,]
    idx <- which(is.na(info_standard$entrezgene_id))
    info_standard <- info_standard[-idx,]
    info_standard[which(info_standard$strand == 1),6] <- "+"
    info_standard[which(info_standard$strand == -1),6] <- "-"
    
    ## Convert to dataframe
    info_dataframe <- data.frame(chr = paste0("chr", info_standard$chromosome_name),
                                 start = info_standard$start_position,
                                 end = info_standard$end_position,
                                 strand = info_standard$strand,
                                 entrez_id = info_standard$entrezgene_id,
                                 hgnc_symbol = c(info_standard$hgnc_symbol))
    
    ## Convert to GRanges
    info_ranges <- makeGRangesFromDataFrame(info_dataframe, keep.extra.columns = TRUE, ignore.strand = FALSE)
    
    ## Return
    return(info_ranges)
    
  } # end generate_all_gene_annotation_as_granges
  
  ####
  get_unique_granges <- function(ranges1, ranges2, min_overlap){
    
    ##
    suppressMessages(library(GenomicRanges))
    
    ##
    overlaps <- findOverlaps(ranges1, ranges2, minoverlap = min_overlap)
    
    ##
    idx_ranges1 <- overlaps@from
    idx_ranges2 <- overlaps@to
    
    ## 
    ranges1_unique <- ranges1[-idx_ranges1]
    ranges2_unique <- ranges2[-idx_ranges2]
    
    ##
    return_ranges <- list(ranges1_unique = ranges1_unique, ranges2_unique = ranges2_unique)
    
    ##
    return(return_ranges)
    
  } # end get_unique_granges
  
  ####
  get_shared_granges <- function(ranges1, ranges2, min_overlap){
    
    ##
    suppressMessages(library(GenomicRanges))
    
    ##
    overlaps <- findOverlaps(ranges1, ranges2, minoverlap = min_overlap)
    
    ##
    idx_ranges1 <- overlaps@from
    idx_ranges2 <- overlaps@to
    
    ## 
    ranges1_shared <- ranges1[idx_ranges1]
    ranges2_shared <- ranges2[idx_ranges2]
    
    ##
    ranges_list <- list(ranges1_shared, ranges2_shared)
    
    ##
    shared_ranged <- get_granges_union(ranges_list)
    
    ##
    return(shared_ranged)
    
  } # end get_shared_granges
  
}
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#############################################################################################################################################################################
#### HOMER ##################################################################################################################################################################
#############################################################################################################################################################################
{
  
  ####
  export_granges_as_HOMER_peaks <- function(ranges, outpath){
    
    ##
    numPeaks <- length(ranges)
    name <- 1:numPeaks
    
    ## retrieve required data for homer peaks
    homerPeaks <- data.frame(name = name,
                             chr = as.character(ranges@seqnames),
                             start = ranges@ranges@start,
                             end = (ranges@ranges@start + ranges@ranges@width),
                             strand = 0)
    
    ## write to file
    write.table(homerPeaks,
                file = outpath,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
  } # export_granges_as_HOMER_peaks
  
  
}
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#############################################################################################################################################################################
#### TOBIAS #################################################################################################################################################################
#############################################################################################################################################################################
{
  
  ####
  export_granges_as_TOBIAS_peaks <- function(ranges, outpath, sample_name){
    
    ##
    chrs <- as.character(ranges@seqnames)
    starts <- ranges@ranges@start
    ends <- ranges@ranges@start + ranges@ranges@width
    names <- rep(sample_name, times = length(chrs))
    
    ##
    df_out <- data.frame(first = chrs, second = starts, third = ends, fourth = names)
    
    ##
    write.table(df_out,
                file = outpath,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    
  } # end export_granges_as_TOBIAS_peaks
  
}
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################





