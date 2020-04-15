
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

####
removeOutlierWidths <- function(inputWidths){
  outlierWidths <- boxplot(inputWidths, plot = FALSE)$out
  outlierIdx <- which(inputWidths %in% outlierWidths)
  inputWidths <- inputWidths[-outlierIdx]
  return(inputWidths)}

####
getSignificantEnrichmentsCompressor <- function(list.enriched, pvalue = 0.05){
  
  ##
  molecular_functions <- as_tibble(list.enriched[["GO_Molecular_Function_2018"]])
  cellular_component <- as_tibble(list.enriched[["GO_Cellular_Component_2018"]])
  biological_process <- as_tibble(list.enriched[["GO_Biological_Process_2018"]])
  genome_browser_PWM <- as_tibble(list.enriched[["Genome_Browser_PWMs"]])
  transfac_jaspar_PWM <- as_tibble(list.enriched[["TRANSFAC_and_JASPAR_PWMs"]])
  encode_tf_chip <- as_tibble(list.enriched[["ENCODE_TF_ChIP-seq_2015"]])
  wikipathways <- as_tibble(list.enriched[["WikiPathways_2019_Human"]])
  kegg <- as_tibble(list.enriched[["KEGG_2019_Human"]])
  tf_lof <- as_tibble(list.enriched[["TF-LOF_Expression_from_GEO"]])
  epigenomics_roadmap <- as_tibble(list.enriched[["Epigenomics_Roadmap_HM_ChIP-seq"]])
  omim_disease <- as_tibble(list.enriched[["OMIM_Disease"]])
  tissue_protein_expression <- as_tibble(list.enriched[["Tissue_Protein_Expression_from_ProteomicsDB"]])
  encode_histone_modifications <- as_tibble(list.enriched[["ENCODE_Histone_Modifications_2015"]])
  
  ####
  sig.molecular_functions <- filter(molecular_functions, Adjusted.P.value <= pvalue)
  sig.cellular_component <- filter(cellular_component, Adjusted.P.value <= pvalue)
  sig.biological_process <- filter(biological_process, Adjusted.P.value <= pvalue)
  sig.genome_browser_PWM  <- filter(genome_browser_PWM, Adjusted.P.value <= pvalue)
  sig.transfac_jaspar_PWM <- filter(transfac_jaspar_PWM, Adjusted.P.value <= pvalue)
  sig.encode_tf_chip <- filter(encode_tf_chip, Adjusted.P.value <= pvalue)
  sig.wikipathways <- filter(wikipathways, Adjusted.P.value <= pvalue)
  sig.kegg <- filter(kegg, Adjusted.P.value <= pvalue)
  sig.tf_lof <- filter(tf_lof, Adjusted.P.value <= pvalue)
  sig.epigenomics_roadmap <- filter(epigenomics_roadmap, Adjusted.P.value <= pvalue)
  sig.omim_disease <- filter(omim_disease, Adjusted.P.value <= pvalue)
  sig.tissue_protein_expression <- filter(tissue_protein_expression, Adjusted.P.value <= pvalue)
  sig.encode_histone_modifications <- filter(encode_histone_modifications, Adjusted.P.value <= pvalue)
  
  #### Return list
  list.return <- list()
  list.return$GO_Molecular_Function_2018 <- sig.molecular_functions
  list.return$GO_Cellular_Component_2018 <- sig.cellular_component
  list.return$GO_Biological_Process_2018 <- sig.biological_process
  list.return$Genome_Browser_PWMs <- sig.genome_browser_PWM
  list.return$TRANSFAC_and_JASPAR_PWMs <- sig.transfac_jaspar_PWM
  list.return$ENCODE_TF_ChIP_seq_2015 <- sig.encode_tf_chip
  list.return$WikiPathways_2019_Human <- sig.wikipathways
  list.return$KEGG_2019_Human <- sig.kegg
  list.return$TF_LOF_Expression_from_GEO <- sig.tf_lof
  list.return$Epigenomics_Roadmap_HM_ChIP_seq <- sig.epigenomics_roadmap
  list.return$OMIM_Disease <- sig.omim_disease
  list.return$Tissue_Protein_Expression_from_ProteomicsDB <- sig.tissue_protein_expression
  list.return$ENCODE_Histone_Modifications_2015 <- sig.encode_histone_modifications
  
  #### Return
  return(list.return)
  
}

#### enrichR compressor
enrichrCompressorPromoters <- function(inputPeaks, pvalue = 0.1){
  
  ##
  dbs.enrichr <- c("GO_Molecular_Function_2018",
                   "GO_Cellular_Component_2018",
                   "GO_Biological_Process_2018",
                   "Genome_Browser_PWMs",
                   "TRANSFAC_and_JASPAR_PWMs",
                   "ENCODE_TF_ChIP-seq_2015",
                   "WikiPathways_2019_Human",
                   "KEGG_2019_Human",
                   "TF-LOF_Expression_from_GEO",
                   "Epigenomics_Roadmap_HM_ChIP-seq",
                   "OMIM_Disease",
                   "Tissue_Protein_Expression_from_ProteomicsDB",
                   "ENCODE_Histone_Modifications_2015")
  
  ##
  temp.peaks.split <- splitGRangesToPromoterAndDistal(inputPeaks)
  temp.peaks.promoters <- temp.peaks.split[["promoterPeaks"]]
  
  ##
  temp.anno <- annotatePeakCompressor(temp.peaks.promoters)
  temp.genes <- unique(temp.anno@anno@elementMetadata@listData[["SYMBOL"]])
  temp.enrichment <- enrichr(temp.genes, dbs.enrichr)
  temp.significant <- getSignificantEnrichmentsCompressor(temp.enrichment, pvalue = pvalue)
  
  ##
  return(temp.significant)
  
}

####
convertToTermsListCompressor <- function(inputTable){
  temp.terms <- c(inputTable[,1])
  temp.terms <- temp.terms[["Term"]]
  return(temp.terms)}

####
removeCommonCompressor <- function(input1, input2){
  
  idx <- which(input1 %in% input2)
  output <- input1[-idx]
  return(output)}

####
importCleanPeaksCompressor <- function(peakPath){
  
  ## Import
  temp.peaks <- readNarrowPeak(peakPath)
  
  ## Cleanup
  temp.peaks <- removeDuplicateGR(cleanGRangesBasic(temp.peaks))
  
  ## Return
  return(temp.peaks)
}

#### makeTSNEvariablityPlot
makeTSNEvariablityPlot <- function(variability, deviation, threshold, perplexity, gene, outpath){
  
  ##
  p0 <- plotVariability(variability, use_plotly = FALSE)
  
  ##
  sample_cor <- getSampleCorrelation(deviation)
  
  ##
  set.seed(seed)
  
  ##
  tsne_results <- deviationsTsne(deviation, threshold = threshold, perplexity = perplexity)
  
  ##
  tsne_plots <- plotDeviationsTsne(deviation, 
                                   tsne_results, 
                                   annotation_name = gene, 
                                   sample_column = "condition", 
                                   shiny = FALSE)
  
  ##
  p1 <- tsne_plots[[1]]
  p2 <- tsne_plots[[2]]
  plot_bottom_row <- plot_grid(p1, p2,ncol = 2)
  p <- plot_grid(p0,
                 plot_bottom_row, 
                 labels = c('A', 'B'),
                 nrow = 2,
                 ncol = 1)
  
  ##
  pdf(file.path(outpath))
  print(p)
  dev.off()	
  
} # end makeTSNEvariablityPlot function
