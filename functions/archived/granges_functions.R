
#### Generate Granges with bin for all genes, promoters, exons, introns
generateBinsPromotersExonsIntronsDistal <- function(promoterExtension = 500){
  
  #### Load libraries ####
  {
    
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    suppressMessages(library(org.Hs.eg.db))
    suppressMessages(library(GenomicRanges))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
  }
  
  #### Promoter annotations ####
  {
    
    ## get the promoters
    promoters <- promoters(txdb, upstream = promoterExtension, downstream = 0, use.names = TRUE)
    
    ## retrieve/update the ENTREZ gene ID and gene symbol for all features, add feature types
    txid <- promoters@elementMetadata@listData[["tx_name"]]
    entrezID <- as.character(mapIds(txdb, keys = txid, column = "GENEID", keytype = "TXNAME"))
    symbolID <- as.character(mapIds(org.Hs.eg.db, keys = entrezID, column = "SYMBOL", keytype = "ENTREZID"))
    ##
    mcols(promoters) <- NULL
    mcols(promoters)$tx_name <- txid
    mcols(promoters)$entrez <- entrezID
    mcols(promoters)$symbol <- symbolID
    mcols(promoters)$feature <- "promoter"
    promoters <- sort(promoters)
    
  }
  
  #### Exon annotations ####
  {
    
    ## get the exons
    exons <- exonsBy(txdb, use.names = TRUE)
    exons <- unlist(exons, use.names = TRUE)
    
    ## retrieve/update the ENTREZ gene ID and gene symbol for all features, add feature types
    txid <- exons@ranges@NAMES
    entrezID <- as.character(mapIds(txdb, keys = txid, column = "GENEID", keytype = "TXNAME"))
    symbolID <- as.character(mapIds(org.Hs.eg.db, keys = entrezID, column = "SYMBOL", keytype = "ENTREZID"))
    ##
    mcols(exons) <- NULL
    mcols(exons)$tx_name <- txid
    mcols(exons)$entrez <- entrezID
    mcols(exons)$symbol <- symbolID
    mcols(exons)$feature <- "exon"
    exons <- sort(exons)
    
  }
  
  #### Intron annotations ####
  {
    
    ## get the introns
    introns <- intronsByTranscript(txdb, use.names = TRUE)
    introns <- unlist(introns, use.names = TRUE)
    
    ## retrieve/update the ENTREZ gene ID and gene symbol for all features, add feature types
    txid <- introns@ranges@NAMES
    entrezID <- as.character(mapIds(txdb, keys = txid, column = "GENEID", keytype = "TXNAME"))
    symbolID <- as.character(mapIds(org.Hs.eg.db, keys = entrezID, column = "SYMBOL", keytype = "ENTREZID"))
    ##
    mcols(introns) <- NULL
    mcols(introns)$tx_name <- txid
    mcols(introns)$entrez <- entrezID
    mcols(introns)$symbol <- symbolID
    mcols(introns)$feature <- "intron"
    introns <- sort(introns)
    
  }
  
  #### Merge promoters, exons, introns ####
  {
    
    ## merge them
    annotations <- c(promoters, exons, introns)
    annotations <- sort(annotations)
    
  }
  
  #### Generate the distal annotations ####
  {
    
    ## now, the the distal regions, generate a GRanges of hg38 and remove anything in the annotations object
    hg38info <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
    hg38names <- as.character(hg38info@seqnames)
    hg38starts <- rep(1, times = 455)
    hg38ends <- hg38info@seqlengths
    hg38df <- data.frame(chromosome = hg38names,
                         start = hg38starts,
                         end = hg38ends)
    ##
    gr.hg38 <- makeGRangesFromDataFrame(hg38df)
    gr.hg38 <- keepStandardChromosomes(gr.hg38, pruning.mode = "coarse")
    gr.hg38 <- dropSeqlevels(gr.hg38, "chrM", pruning.mode = "coarse")
    gr.hg38 <- sort(gr.hg38)
    
    ## now get the set difference
    distal <- setdiff(gr.hg38, annotations, ignore.strand = TRUE)
    mcols(distal) <- NULL
    mcols(distal)$tx_name <- "distal"
    mcols(distal)$entrez <- "distal"
    mcols(distal)$symbol <- "distal"
    mcols(distal)$feature <- "distal"
    numDistal <- length(distal)
    distal@ranges@NAMES <- as.character(rep("distal", times = numDistal))
    
  }
  
  #### Merge to get the final complete annotation object ####
  {
    
    ## merge for the final annotation object
    final <- c(annotations, distal)
    final <- sortSeqlevels(final)
    final <- sort(final, ignore.strand = TRUE)
    
  }

} # end generateBinsPromotersExonsIntronsDistal

####
generateAllGeneAnnotationsAsGRanges <- function(){
  
  ##
  library(biomaRt)
  library(GenomicRanges)
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
  
}

####
generateStandardChromosomeNames <- function(){
  names <- paste0("chr", c(1:22, "X", "Y"))
  return(names)
}

#### Split a GRanges to promoter and distal subset, distance defines region around TSS that is considered a promoter ####
splitGRangesToPromoterAndDistal <- function(ranges, distance = 500){
  
  ##
  suppressMessages(library(GenomicRanges))
  
  ##
  promoterWindows <- getGenePromoterWindows(window = distance)
  
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
  
} # end splitGRangesToPromoterAndDistal

#### Remove duplicate GRanges ####
removeDuplicateGRanges <- function(ranges){
  
  suppressMessages(library(GenomicRanges))
  deduplicated <- unique(ranges)
  return(deduplicated)
  
} # end removeDuplicateGR

#### Remove duplicate GRanges ####
removeOverlappingGRanges <- function(subjectRanges, queryRanges){
  
  ####
  suppressMessages(library(GenomicRanges))
  
  tempRemove <- subsetByOverlaps(subjectRanges, queryRanges)
  tempKeep <- GenomicRanges::setdiff(subjectRanges, tempRemove, ignore.strand = TRUE)
  return(tempKeep)
  
} # end removeDuplicateGR

#### Merge two or more GRanges objects ####
getGRangesUnion <- function(grList = list()){
  
  suppressMessages(library(GenomicRanges))
  temp.merged <- grList[[1]]
  numItems <- length(grList)
  
  ##
  for (a in 2:numItems){
    temp.gr <- grList[[a]]
    temp.merged <- GenomicRanges::union(temp.merged, temp.gr)
  }
  
  ##
  return(temp.merged)
  
} # end getGRangesUnion

#### Make genomic windows as a GRanges object ####
generateGenomicWindows <- function(binSize = 100, clean = TRUE){
  
  ##
  cat("Loading libraries for generateGenomicWindows function \n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  
  ## Important to force use of hg38 in case hg19 loaded also
  genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  
  ##
  bins = tileGenome(seqinfo(genome),
                    tilewidth = binSize,
                    cut.last.tile.in.chrom = TRUE)
  
  ## Clean if true
  if (clean == TRUE){bins <- cleanGRanges(bins)}
  
  ##
  return(bins)
  
}






