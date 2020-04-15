
#### GVIZ plotting validation to check consistency between BAM, peak, summit, and bigwig files
validate.gviz_plot_bam_peak_summit_bigwig <- function(bamPath, peakPath, bigwigPath, plotBin, sampleName, verbose = TRUE){
  
  ## Load libraries ##
  {
    suppressMessages(library(Gviz))
    suppressMessages(library(GenomicRanges))
  }
  
  ## Get the plotting details ##
  {
    
    ##
    genome = "hg38"
    
    ##
    currentChr <- as.character(plotBin@seqnames)
    currentStart <- plotBin@ranges@start
    currentEnd <- currentStart + plotBin@ranges@width
    
    ##
    if(verbose){cat("Current chromosome:", currentChr, "\n")
      cat("Current start:", currentStart, "\n")
      cat("Current end:", currentEnd, "\n")
    }}
  
  ## Generate the BAM track ##
  {
    
    if(verbose){cat("Generating bam track \n")}
    ##
    bamTrack <- AlignmentsTrack(range = bamPath,
                                genome = genome,
                                chromosome = currentChr,
                                start = currentStart,
                                end = currentEnd,
                                strand = "*",
                                #cex = at.cex,
                                type = "coverage",
                                showTitle = TRUE,
                                name = "BAM",                 # title of the track
                                #ylim = at.ylim,
                                #fontcolor.title = titleColor, # color of the title
                                col.coverage = "black", # the outline color of the coverage lines
                                fill.coverage = "black") # the fill color of the coverage body
    
  }
  
  ## Generate the peak track ##
  {
    
    if(verbose){cat("Generating peak track \n")}
    ##
    peaks <- importBED(peakPath)
    ##
    peakTrack <- DataTrack(range = peaks,
                           genome = genome,
                           #name = "PEAKS",
                           type = "histogram",                             ## Type of data track (histogram here)
                           #ylim = dt.ylim,                             ## Range of the y-axis
                           #lwd = dt.lwd,                               ## Default line width scaling
                           fontcolor.title = "black",                  ## Font color of the title
                           col.histogram = "black",             ## Histogram color
                           fill = "black"                           ## Fill color
                           #lwd.border.title = dt.lwd.border.title,
                           #lineheight = dt.lineheight,
                           #fontsize = dt.fontsize,
                           #frame = dt.frame,
                           #cex.title = dt.cex.title,
                           #rotation.title = dt.rotation.title,
                           #lwd.title = dt.lwd.title,
                           #fontface.title = dt.fontface.title,
                           #col.frame = dt.col.frame
                           )
    
  }
  
  ## Generate the bigwig track ##
  {
    
    if(verbose){cat("Generating bigwig track \n")}
    ##
    bigwigTrack <- DataTrack(range = bigwigPath,
                           genome = genome,
                           name = "BIGWIG",
                           type = "histogram",                             ## Type of data track (histogram here)
                           #ylim = dt.ylim,                             ## Range of the y-axis
                           #lwd = dt.lwd,                               ## Default line width scaling
                           fontcolor.title = "black",                  ## Font color of the title
                           col.histogram = "black",             ## Histogram color
                           fill = "black"                           ## Fill color
                           #lwd.border.title = dt.lwd.border.title,
                           #lineheight = dt.lineheight,
                           #fontsize = dt.fontsize,
                           #frame = dt.frame,
                           #cex.title = dt.cex.title,
                           #rotation.title = dt.rotation.title,
                           #lwd.title = dt.lwd.title,
                           #fontface.title = dt.fontface.title,
                           #col.frame = dt.col.frame
    ) 
    
  }
  
  #### Make the bin track ####
  {
    track.BIN <- makeBigwigTrack(plotBin, " ", "red", "red")
  }
  
  #### Generate the plot ####
  {
    
    ##
    plotTracks(list(bamTrack, peakTrack, bigwigTrack),
               from = currentStart - 1000,
               to = currentEnd + 1000,
               sizes = c(1,1,1),
               main = sampleName,
               cex.main = 1)
    
  }
  
  
  
}

#### Compare BAM reads to FASTQ reads
validate.BAM_FASTQ_barcode_match <- function(inputBED, outputPath){
  
}

#### Validate all objects in a list are unique ####
validate.check_all_objects_unique <- function(inputList, verbose = FALSE){
  
  ##
  message(">>> Performing all objects unique validation <<<")
  
  ##
  unique = TRUE
  
  ##
  numObjects <- length(inputList)
  if(verbose){cat("Detected", numObjects, "input objects \n")}
  
  ##
  for (a in 1:numObjects){
    
    ## Get the subject object
    subject <- inputList[[a]]
    
    ## Get the query objects
    queryList <- inputList[-a]
    numQuery <- length(queryList)
    
    ##
    for (b in 1:numQuery){
      
      ##
      com <- paste0("query <- queryList[[", b, "]]")
      eval(parse(text = com))
      
      ##
      if(identical(subject, query)){identical = TRUE}
      
    }
    
  }
  
  ##
  if(!unique){message(">>> Found duplicated objects!!! <<<")}
  if(unique){message(">>> No duplicate objects found <<<")}
  
  ## return
  return(unique)
  
} # end validate.check_all_objects_unique

#### Validate by plotting peaks in bins as heatmap in linear genome representation ####

