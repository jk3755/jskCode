
#### Set the track and display parameters ####
setParameters <- function(
  #### Locus ####
  genome,
  chromosome,
  start,
  end,
  strand,
  #### Alignment track ####
  at.showTitle,
  at.cex,
  at.ylim,
  at.type,
  #### Data track ####
  dt.type,
  dt.ylim,
  dt.lwd,
  dt.titlecolor,
  dt.lwd.border.title,
  dt.lineheight,
  dt.fontsize,
  dt.frame,
  dt.cex.title,
  dt.rotation.title,
  dt.lwd.title,
  dt.fontface.title,
  dt.col.frame,
  #### Gene track ####
  gt.showID,
  gt.geneSymbol,
  gt.symbol,
  #### Plot ####
  plot.sizes,
  plot.upstream,
  plot.downstream,
  plot.title){
  
  #### Locus ####
  {
    genome <<- genome
    chromosome <<- chromosome
    start <<- start
    end <<- end
    strand <<- strand
  }
  
  #### Alignment track options ####
  {
    at.showTitle <<- at.showTitle
    at.cex <<- at.cex
    at.ylim <<- at.ylim
    at.type <<- at.type
  }
  
  #### Data track options ####
  {
    dt.type <<- dt.type
    dt.ylim <<- dt.ylim
    dt.lwd <<- dt.lwd
    dt.titlecolor <<- dt.titlecolor
    dt.lwd.border.title <<- dt.lwd.border.title
    dt.lineheight <<- dt.lineheight
    dt.fontsize <<- dt.fontsize
    dt.frame <<- dt.frame
    dt.cex.title <<- dt.cex.title
    dt.rotation.title <<- dt.rotation.title
    dt.lwd.title <<- dt.lwd.title
    dt.fontface.title <<- dt.fontface.title
    dt.col.frame <<- dt.col.frame
  }
  
  #### Plotting options ####
  {
    plot.sizes <<- plot.sizes                     ## Numeric vector, the relative size of each track
    plot.from <<- start - plot.upstream           ## Starting position of the plot
    plot.to <<- end + plot.downstream             ## End position of the plot
    plot.title <<- plot.title                     ## Title of the plot
  }
  
} # end setParameters function

#### Generate BIGWIG tracks (datatrack) ####
makeBigwigTrack <- function(bigwig, trackName, histogramColor, fillColor){
  
  ## displayPars() allows you to display all parameters for a given track
  
  ##
  cat("Generating bigwig track with the following parameters...", "\n")
  cat("Track name:", trackName, "\n")
  cat("Histogram color:", histogramColor, "\n")
  
  ## Initialize the track
  cat("Initializing DataTrack", "\n")
  track.bigwig <- DataTrack(range = bigwig,
                            genome = genome,
                            name = trackName,
                            type = dt.type,                             ## Type of data track (histogram here)
                            ylim = dt.ylim,                             ## Range of the y-axis
                            lwd = dt.lwd,                               ## Default line width scaling
                            fontcolor.title = "black",                  ## Font color of the title
                            col.histogram = histogramColor,             ## Histogram color
                            fill = fillColor,                           ## Fill color
                            lwd.border.title = dt.lwd.border.title,
                            lineheight = dt.lineheight,
                            fontsize = dt.fontsize,
                            frame = dt.frame,
                            cex.title = dt.cex.title,
                            rotation.title = dt.rotation.title,
                            lwd.title = dt.lwd.title,
                            fontface.title = dt.fontface.title,
                            col.frame = dt.col.frame)         
  
  ## Return the track
  cat("Returning DataTrack", "\n")
  return(track.bigwig)
  
} # end makeBigwigTrack function

#### Generate BIGWIG plots ####
generatePlotsBigwigAveraged <- function(chr, start, end, plot.upstream, plot.downstream, name, outputDir, plotWidth, plotHeight, plotDPI){
  
  ## Report ##
  {
    cat("Generating Bigwig plot", "\n")
    cat("Current chromosome:", chr, "\n")
    cat("Current start:", start, "\n")
    cat("Current end:", end, "\n")
    cat("Current plot name:", name, "\n")
    cat("Output directory for figures:", outputDir, "\n")
    cat("Plot width:", plotWidth, "\n")
    cat("Plot height:", plotHeight, "\n")
    cat("Plot DPI:", plotDPI, "\n")
  }
  
  ## Set parameters ##
  {
    cat("Setting track parameters", "\n")
    setParameters(
      #### Locus ####
      genome = "hg38",
      chromosome = chr,
      start = start,
      end = end,
      strand = "*",
      #### Alignment track ####
      at.showTitle = TRUE,
      at.cex = 0.5,                
      at.ylim = c(0,5),
      at.type = "coverage",
      #### Data track ####
      dt.type = "histogram",            ## Type of track for the bigwig
      dt.ylim = c(0,4),                 ## Yaxis range for bigwig
      dt.lwd = 1,
      dt.lineheight = 0.5,
      dt.fontsize = 9,
      dt.frame = TRUE,
      dt.cex.title = 0.75,
      dt.titlecolor = "black",
      dt.lwd.border.title = 1,
      dt.rotation.title = 90,
      dt.lwd.title = 1,
      dt.fontface.title = 1,
      dt.col.frame = "black",
      #### Plotting ####
      plot.sizes = c(
        0.5,                            # Ideogram track
        1,                              # PROG track
        1,                              # BAZ2B track
        1,                              # LUF track
        1,                              # HSC track
        1,                              # MPP track
        1,                              # LMP track
        1,                              # CMP track
        1,                              # GMP track
        1,                              # MEP track
        1,                              # MONO track
        1,                              # ERY track
        1,                              # CLP track
        1,                              # CD4 track
        1,                              # CD8 track
        1,                              # B track
        1,                              # NK track
        1.5),                             # Gene annotation track
      plot.upstream = plot.upstream,
      plot.downstream = plot.downstream,
      plot.title = name
    )
  }
  
  ####
  cat("Initializing tracks", "\n")
  track.ideogram <- IdeogramTrack(genome = genome, chromosome = chromosome, start = start, end = end)
  
  #### The gene annotations track with gene symbols ####
  {
    ucscGenes <- UcscTrack(genome = genome,
                           chromosome = chromosome,
                           table = "ncbiRefSeq",
                           track = 'NCBI RefSeq',
                           trackType = "GeneRegionTrack",
                           rstarts = "exonStarts",
                           rends = "exonEnds",
                           gene = "name",
                           symbol = 'name',
                           transcript = "name",
                           strand = "strand",
                           stacking = 'hide',
                           name = "GENE",
                           showID = TRUE,
                           geneSymbol = TRUE,
                           showTitle = TRUE,
                           fontcolor.title = "black",        ## Color of the track title font
                           lineheight = 0.5,                 ## Font height of all text
                           rotation.group = 0,               ## Degree of text rotation for group labels
                           just.group = "left",              ## Justification of group labels
                           fontsize.group = 10,              ## Font size for group level annotation
                           frame = TRUE,
                           col.frame = "black",
                           cex.axis = 0.5,
                           fontsize = 12,
                           fontcolor = "black",
                           col.symbol = "black",
                           fontface.group = 2,
                           fontcolor.group = "black",
                           fontcolor.item = "black"
    )
  }
  
  #### Add the gene symbols to the track ####
  z <- ranges(ucscGenes)
  mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
  ucscGenes2 <- ucscGenes
  ranges(ucscGenes2) <- z
  track.UCSCgenes.symbols <- ucscGenes2
  
  #### Initiate the Bigwig tracks ####
  {
    cat("Initializing Bigwig tracks", "\n")
    
    ##
    # Filenames validated []
    #
    track.PROG <- makeBigwigTrack(avg.bw.PROG, "\n\nPROG\n\n(n = 2)", "black")
    #
    track.BAZ2B <- makeBigwigTrack(avg.bw.BAZ2B, "\n\nBAZ2B\n\n(n = 2)", "forestgreen")
    #
    track.LUF <- makeBigwigTrack(avg.bw.LUF, "\n\nLUF\n\n(n = 2)", "greenyellow")
    #
    track.HSC <- makeBigwigTrack(avg.bw.HSC, "\n\nHSC\n\n(n = 7)", "darkgreen")
    #
    track.MPP <- makeBigwigTrack(avg.bw.MPP, "\n\nMPP\n\n(n = 5)", "green3")
    #
    track.LMP <- makeBigwigTrack(avg.bw.LMP, "\n\nLMP\n\n(n = 3)", "cyan2")
    #
    track.CMP <- makeBigwigTrack(avg.bw.CMP, "\n\nCMP\n\n(n = 6)", "peachpuff")
    #
    track.GMP <- makeBigwigTrack(avg.bw.GMP, "\n\nGMP\n\n(n = 7)", "goldenrod1")
    #
    track.MEP <- makeBigwigTrack(avg.bw.MEP, "\n\nMEP\n\n(n = 6)", "red1")
    #
    track.MONO <- makeBigwigTrack(avg.bw.MONO, "\n\nMONO\n\n(n = 6)", "orangered1")
    #
    track.ERY <- makeBigwigTrack(avg.bw.ERY, "\n\nERY\n\n(n = 6)", "indianred4")
    #
    track.CLP <- makeBigwigTrack(avg.bw.CLP, "\n\nCLP\n\n(n = 4)", "cadetblue1")
    #
    track.CD4 <- makeBigwigTrack(avg.bw.CD4, "\n\nCD4\n\n(n = 4)", "cornflowerblue")
    #
    track.CD8 <- makeBigwigTrack(avg.bw.CD8, "\n\nCD8\n\n(n = 4)", "blue3")
    #
    track.B <- makeBigwigTrack(avg.bw.B, "\n\nB\n\n(n = 3)", "darkorchid1")
    #
    track.NK <- makeBigwigTrack(avg.bw.NK, "\n\nNK\n\n(n = 5)", "darkorchid4")
    
  }
  
  #### Generate the plot ###
  cat("Generating plot", "\n")
  figureOutpath <- file.path(outputDir, paste0(name, ".png"))
  cat("Figure output filepath:", figureOutpath, "\n")
  png(filename = figureOutpath, width = plotWidth, height = plotHeight, res = plotDPI)
  
  #### Generate the track plot ####
  {
    plotTracks(list(track.ideogram,
                    track.PROG,
                    track.BAZ2B,
                    track.LUF,
                    track.HSC,
                    track.MPP,
                    track.LMP,
                    track.CMP,
                    track.GMP,
                    track.MEP,
                    track.MONO,
                    track.ERY,
                    track.CLP,
                    track.CD4,
                    track.CD8,
                    track.B,
                    track.NK,
                    track.UCSCgenes.symbols),
               from = plot.from,
               to = plot.to,
               sizes = plot.sizes,
               main = plot.title,
               cex.main = 1,
               transcriptAnnotation = "symbol",
               title.width = 0.5)
  }
  
  ##
  dev.off()
  
} # end generatePlotsBigwig function


#### Generate BIGWIG plots ####
generatePlotsPeaksAveraged <- function(chr, start, end, plot.upstream, plot.downstream, name, outputDir, plotWidth, plotHeight, plotDPI){
  
  ## Report ##
  {
    cat("Generating Bigwig plot", "\n")
    cat("Current chromosome:", chr, "\n")
    cat("Current start:", start, "\n")
    cat("Current end:", end, "\n")
    cat("Current plot name:", name, "\n")
    cat("Output directory for figures:", outputDir, "\n")
    cat("Plot width:", plotWidth, "\n")
    cat("Plot height:", plotHeight, "\n")
    cat("Plot DPI:", plotDPI, "\n")
  }
  
  ## Set parameters ##
  {
    cat("Setting track parameters", "\n")
    setParameters(
      #### Locus ####
      genome = "hg38",
      chromosome = chr,
      start = start,
      end = end,
      strand = "*",
      #### Alignment track ####
      at.showTitle = TRUE,
      at.cex = 0.5,                
      at.ylim = c(0,5),
      at.type = "coverage",
      #### Data track ####
      dt.type = "histogram",            ## Type of track for the bigwig
      dt.ylim = c(0,100),                 ## Yaxis range for bigwig
      dt.lwd = 1,
      dt.lineheight = 0.5,
      dt.fontsize = 9,
      dt.frame = TRUE,
      dt.cex.title = 0.75,
      dt.titlecolor = "black",
      dt.lwd.border.title = 1,
      dt.rotation.title = 90,
      dt.lwd.title = 1,
      dt.fontface.title = 1,
      dt.col.frame = "black",
      #### Plotting ####
      plot.sizes = c(
        0.5,                            # Ideogram track
        0.25,                           # BIN track
        1,                              # PROG track
        1,                              # BAZ2B track
        1,                              # LUF track
        1,                              # HSC track
        1,                              # MPP track
        1,                              # LMP track
        1,                              # CMP track
        1,                              # GMP track
        1,                              # MEP track
        1,                              # MONO track
        1,                              # ERY track
        1,                              # CLP track
        1,                              # CD4 track
        1,                              # CD8 track
        1,                              # B track
        1,                              # NK track
        1.5),                             # Gene annotation track
      plot.upstream = plot.upstream,
      plot.downstream = plot.downstream,
      plot.title = name
    )
  }
  
  ####
  cat("Initializing tracks", "\n")
  track.ideogram <- IdeogramTrack(genome = genome, chromosome = chromosome, start = start, end = end)
  
  #### The gene annotations track with gene symbols ####
  {
    ucscGenes <- UcscTrack(genome = genome,
                           chromosome = chromosome,
                           table = "ncbiRefSeq",
                           track = 'NCBI RefSeq',
                           trackType = "GeneRegionTrack",
                           rstarts = "exonStarts",
                           rends = "exonEnds",
                           gene = "name",
                           symbol = 'name',
                           transcript = "name",
                           strand = "strand",
                           stacking = 'hide',
                           name = "GENE",
                           showID = TRUE,
                           geneSymbol = TRUE,
                           showTitle = TRUE,
                           fontcolor.title = "black",        ## Color of the track title font
                           lineheight = 0.5,                 ## Font height of all text
                           rotation.group = 0,               ## Degree of text rotation for group labels
                           just.group = "left",              ## Justification of group labels
                           fontsize.group = 10,              ## Font size for group level annotation
                           frame = TRUE,
                           col.frame = "black",
                           cex.axis = 0.5,
                           fontsize = 12,
                           fontcolor = "black",
                           col.symbol = "black",
                           fontface.group = 2,
                           fontcolor.group = "black",
                           fontcolor.item = "black"
    )
  }
  
  #### Add the gene symbols to the track ####
  z <- ranges(ucscGenes)
  mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
  ucscGenes2 <- ucscGenes
  ranges(ucscGenes2) <- z
  track.UCSCgenes.symbols <- ucscGenes2
  
  #### Initiate the Bigwig tracks ####
  {
    cat("Initializing Bigwig tracks", "\n")
    
    ##
    # Filenames validated []
    #
    track.PROG <- makeBigwigTrack(peaks.pooled.PROG, "\n\nPROG\n\n(n = 2)", "black", "black")
    #
    track.BAZ2B <- makeBigwigTrack(peaks.pooled.BAZ2B, "\n\nBAZ2B\n\n(n = 2)", "forestgreen", "black")
    #
    track.LUF <- makeBigwigTrack(peaks.pooled.LUF, "\n\nLUF\n\n(n = 2)", "greenyellow", "black")
    #
    track.HSC <- makeBigwigTrack(peaks.pooled.HSC, "\n\nHSC\n\n(n = 7)", "darkgreen", "black")
    #
    track.MPP <- makeBigwigTrack(peaks.pooled.MPP, "\n\nMPP\n\n(n = 5)", "green3", "black")
    #
    track.LMP <- makeBigwigTrack(peaks.pooled.LMP, "\n\nLMP\n\n(n = 3)", "cyan2", "black")
    #
    track.CMP <- makeBigwigTrack(peaks.pooled.CMP, "\n\nCMP\n\n(n = 6)", "peachpuff", "black")
    #
    track.GMP <- makeBigwigTrack(peaks.pooled.GMP, "\n\nGMP\n\n(n = 7)", "goldenrod1", "black")
    #
    track.MEP <- makeBigwigTrack(peaks.pooled.MEP, "\n\nMEP\n\n(n = 6)", "red1", "black")
    #
    track.MONO <- makeBigwigTrack(peaks.pooled.MONO, "\n\nMONO\n\n(n = 6)", "orangered1", "black")
    #
    track.ERY <- makeBigwigTrack(peaks.pooled.ERY, "\n\nERY\n\n(n = 6)", "indianred4", "black")
    #
    track.CLP <- makeBigwigTrack(peaks.pooled.CLP, "\n\nCLP\n\n(n = 4)", "cadetblue1", "black")
    #
    track.CD4 <- makeBigwigTrack(peaks.pooled.CD4, "\n\nCD4\n\n(n = 4)", "cornflowerblue", "black")
    #
    track.CD8 <- makeBigwigTrack(peaks.pooled.CD8, "\n\nCD8\n\n(n = 4)", "blue3", "black")
    #
    track.B <- makeBigwigTrack(peaks.pooled.B, "\n\nB\n\n(n = 3)", "darkorchid1", "black")
    #
    track.NK <- makeBigwigTrack(peaks.pooled.NK, "\n\nNK\n\n(n = 5)", "darkorchid4", "black")
    
  }
  
  #### Make the bin track ####
  {
    track.BIN <- makeBigwigTrack(plotBins, " ", "red", "red")
  }
  
  #### Generate the plot ###
  cat("Generating plot", "\n")
  figureOutpath <- file.path(outputDir, paste0(name, ".png"))
  cat("Figure output filepath:", figureOutpath, "\n")
  png(filename = figureOutpath, width = plotWidth, height = plotHeight, res = plotDPI)
  
  #### Generate the track plot ####
  {
    plotTracks(list(track.ideogram,
                    track.BIN,
                    track.PROG,
                    track.BAZ2B,
                    track.LUF,
                    track.HSC,
                    track.MPP,
                    track.LMP,
                    track.CMP,
                    track.GMP,
                    track.MEP,
                    track.MONO,
                    track.ERY,
                    track.CLP,
                    track.CD4,
                    track.CD8,
                    track.B,
                    track.NK,
                    track.UCSCgenes.symbols),
               from = plot.from,
               to = plot.to,
               sizes = plot.sizes,
               main = plot.title,
               cex.main = 1,
               transcriptAnnotation = "symbol",
               title.width = 0.5)
  }
  
  ##
  dev.off()
  
} # end generatePlotsBigwig function

####
generatePlotsBWAveraged <- function(chr, start, end, plot.upstream, plot.downstream, name, outputDir, plotWidth, plotHeight, plotDPI){
  
  ## Report ##
  {
    cat("Generating Bigwig plot", "\n")
    cat("Current chromosome:", chr, "\n")
    cat("Current start:", start, "\n")
    cat("Current end:", end, "\n")
    cat("Current plot name:", name, "\n")
    cat("Output directory for figures:", outputDir, "\n")
    cat("Plot width:", plotWidth, "\n")
    cat("Plot height:", plotHeight, "\n")
    cat("Plot DPI:", plotDPI, "\n")
  }
  
  ## Set parameters ##
  {
    cat("Setting track parameters", "\n")
    setParameters(
      #### Locus ####
      genome = "hg38",
      chromosome = chr,
      start = start,
      end = end,
      strand = "*",
      #### Alignment track ####
      at.showTitle = TRUE,
      at.cex = 0.5,                
      at.ylim = c(0,5),
      at.type = "coverage",
      #### Data track ####
      dt.type = "histogram",            ## Type of track for the bigwig
      dt.ylim = c(0,5),                 ## Yaxis range for bigwig
      dt.lwd = 1,
      dt.lineheight = 0.5,
      dt.fontsize = 9,
      dt.frame = TRUE,
      dt.cex.title = 0.75,
      dt.titlecolor = "black",
      dt.lwd.border.title = 1,
      dt.rotation.title = 90,
      dt.lwd.title = 1,
      dt.fontface.title = 1,
      dt.col.frame = "black",
      #### Plotting ####
      plot.sizes = c(
        0.5,                            # Ideogram track
        0.25,                           # BIN track
        1,                              # PROG track
        1,                              # BAZ2B track
        1,                              # LUF track
        1,                              # HSC track
        1,                              # MPP track
        1,                              # LMP track
        1,                              # CMP track
        1,                              # GMP track
        1,                              # MEP track
        1,                              # MONO track
        1,                              # ERY track
        1,                              # CLP track
        1,                              # CD4 track
        1,                              # CD8 track
        1,                              # B track
        1,                              # NK track
        1.5),                             # Gene annotation track
      plot.upstream = plot.upstream,
      plot.downstream = plot.downstream,
      plot.title = name
    )
  }
  
  ####
  cat("Initializing tracks", "\n")
  track.ideogram <- IdeogramTrack(genome = genome, chromosome = chromosome, start = start, end = end)
  
  #### The gene annotations track with gene symbols ####
  {
    ucscGenes <- UcscTrack(genome = genome,
                           chromosome = chromosome,
                           table = "ncbiRefSeq",
                           track = 'NCBI RefSeq',
                           trackType = "GeneRegionTrack",
                           rstarts = "exonStarts",
                           rends = "exonEnds",
                           gene = "name",
                           symbol = 'name',
                           transcript = "name",
                           strand = "strand",
                           stacking = 'hide',
                           name = "GENE",
                           showID = TRUE,
                           geneSymbol = TRUE,
                           showTitle = TRUE,
                           fontcolor.title = "black",        ## Color of the track title font
                           lineheight = 0.5,                 ## Font height of all text
                           rotation.group = 0,               ## Degree of text rotation for group labels
                           just.group = "left",              ## Justification of group labels
                           fontsize.group = 10,              ## Font size for group level annotation
                           frame = TRUE,
                           col.frame = "black",
                           cex.axis = 0.5,
                           fontsize = 12,
                           fontcolor = "black",
                           col.symbol = "black",
                           fontface.group = 2,
                           fontcolor.group = "black",
                           fontcolor.item = "black"
    )
  }
  
  #### Add the gene symbols to the track ####
  z <- ranges(ucscGenes)
  mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
  ucscGenes2 <- ucscGenes
  ranges(ucscGenes2) <- z
  track.UCSCgenes.symbols <- ucscGenes2
  
  #### Initiate the Bigwig tracks ####
  {
    cat("Initializing Bigwig tracks", "\n")
    
    ##
    # Filenames validated []
    #
    track.PROG <- makeBigwigTrack(avg.bw.PROG, "\n\nPROG\n\n(n = 2)", "black", "black")
    #
    track.BAZ2B <- makeBigwigTrack(avg.bw.BAZ2B, "\n\nBAZ2B\n\n(n = 2)", "forestgreen", "black")
    #
    track.LUF <- makeBigwigTrack(avg.bw.LUF, "\n\nLUF\n\n(n = 2)", "greenyellow", "black")
    #
    track.HSC <- makeBigwigTrack(avg.bw.HSC, "\n\nHSC\n\n(n = 7)", "darkgreen", "black")
    #
    track.MPP <- makeBigwigTrack(avg.bw.MPP, "\n\nMPP\n\n(n = 5)", "green3", "black")
    #
    track.LMP <- makeBigwigTrack(avg.bw.LMP, "\n\nLMP\n\n(n = 3)", "cyan2", "black")
    #
    track.CMP <- makeBigwigTrack(avg.bw.CMP, "\n\nCMP\n\n(n = 6)", "peachpuff", "black")
    #
    track.GMP <- makeBigwigTrack(avg.bw.GMP, "\n\nGMP\n\n(n = 7)", "goldenrod1", "black")
    #
    track.MEP <- makeBigwigTrack(avg.bw.MEP, "\n\nMEP\n\n(n = 6)", "red1", "black")
    #
    track.MONO <- makeBigwigTrack(avg.bw.MONO, "\n\nMONO\n\n(n = 6)", "orangered1", "black")
    #
    track.ERY <- makeBigwigTrack(avg.bw.ERY, "\n\nERY\n\n(n = 6)", "indianred4", "black")
    #
    track.CLP <- makeBigwigTrack(avg.bw.CLP, "\n\nCLP\n\n(n = 4)", "cadetblue1", "black")
    #
    track.CD4 <- makeBigwigTrack(avg.bw.CD4, "\n\nCD4\n\n(n = 4)", "cornflowerblue", "black")
    #
    track.CD8 <- makeBigwigTrack(avg.bw.CD8, "\n\nCD8\n\n(n = 4)", "blue3", "black")
    #
    track.B <- makeBigwigTrack(avg.bw.B, "\n\nB\n\n(n = 3)", "darkorchid1", "black")
    #
    track.NK <- makeBigwigTrack(avg.bw.NK, "\n\nNK\n\n(n = 5)", "darkorchid4", "black")
    
  }
  
  #### Make the bin track ####
  {
    track.BIN <- makeBigwigTrack(plotBins, " ", "red", "red")
  }
  
  #### Generate the plot ###
  cat("Generating plot", "\n")
  figureOutpath <- file.path(outputDir, paste0(name, ".png"))
  cat("Figure output filepath:", figureOutpath, "\n")
  png(filename = figureOutpath, width = plotWidth, height = plotHeight, res = plotDPI)
  
  #### Generate the track plot ####
  {
    plotTracks(list(track.ideogram,
                    track.BIN,
                    track.PROG,
                    track.BAZ2B,
                    track.LUF,
                    track.HSC,
                    track.MPP,
                    track.LMP,
                    track.CMP,
                    track.GMP,
                    track.MEP,
                    track.MONO,
                    track.ERY,
                    track.CLP,
                    track.CD4,
                    track.CD8,
                    track.B,
                    track.NK,
                    track.UCSCgenes.symbols),
               from = plot.from,
               to = plot.to,
               sizes = plot.sizes,
               main = plot.title,
               cex.main = 1,
               transcriptAnnotation = "symbol",
               title.width = 0.5)
  }
  
  ##
  dev.off()
  
} # end generatePlotsBigwig function


#### Generate BIGWIG plots ####
generatePlotsBigwig <- function(chrList, startList, endList, nameList, outputDir, plotWidth, plotHeight, plotDPI, ylim = c(0,10)){
  
  ## Report
  numItems <- length(chrList)
  cat("Generating BIGWIG plots", "\n")
  cat("Number of items:", numItems, "\n")
  cat("Output directory for figures:", outputDir, "\n")
  cat("Plot width:", plotWidth, "\n")
  cat("Plot height:", plotHeight, "\n")
  cat("Plot DPI:", plotDPI, "\n")
  
  ## Loop
  for (a in 1:numItems){
    
    tryCatch({
      
      ## Get current data
      chr <- chrList[a]
      start <- startList[a]
      end <- endList[a]
      name <- nameList[a]
      
      ##
      cat("Generating plot number:", a, "\n")
      cat("Current chromosome:", chr, "\n")
      cat("Current start:", start, "\n")
      cat("Current end:", end, "\n")
      cat("Current gene:", name, "\n")
      
      ## Set paramters for track initiation
      cat("Setting track parameters", "\n")
      setParameters("hg38", chr, start, end, "*",
                    TRUE, 0.5, c(0,5), "coverage",
                    "histogram", ylim, 1, "black")
      
      ##
      cat("Initializing tracks", "\n")
      track.ideogram <- IdeogramTrack(genome = genome, chromosome = chromosome, start = start, end = end)
      track.axis <- GenomeAxisTrack(genome = genome, chromosome = chromosome, start = start, end = end)
      
      ## The gene annotations track with gene symbols
      ucscGenes <- UcscTrack(genome=genome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                             chromosome=chromosome, rstarts = "exonStarts", rends = "exonEnds",
                             gene = "name", symbol = 'name', transcript = "name",
                             strand = "strand", stacking = 'hide', showID = T, geneSymbol = T,
                             showTitle = T, name = "GENE", fontcolor.title = "white",
                             lineheight=1, rotation.group=90, just.group="left", fontsize.group=12)
      z <- ranges(ucscGenes)
      mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
      ucscGenes2 <- ucscGenes
      ranges(ucscGenes2) <- z
      track.UCSCgenes.symbols <- ucscGenes2
      
      #### Initiate the BIGWIG tracks ####
      cat("Initializing BIGWIG tracks", "\n")
      track.bigwig.prog_A <- makeBigwigTrack(bigwig.PROG_A, "PROG_A", "blue")
      track.bigwig.prog_B <- makeBigwigTrack(bigwig.PROG_B, "PRoG_B", "blue")
      track.bigwig.luf_A <- makeBigwigTrack(bigwig.LUF_A, "LUF_A", "yellow")
      track.bigwig.luf_B <- makeBigwigTrack(bigwig.LUF_B, "LUF_B", "yellow")
      track.bigwig.baz2b_A <- makeBigwigTrack(bigwig.BAZ2B_A, "BAZ2B_A", "red")
      track.bigwig.baz2b_B <- makeBigwigTrack(bigwig.BAZ2B_B, "BAZ2B_B", "red")
      
      ## Set plotting parameters
      cat("Setting plotting parameters", "\n")
      setPlotParameters(c(0.5,0.5,1,1,1,1,1,1,2), 5000, 5000, name)
      
      ## Generate the plot
      cat("Generating plot", "\n")
      figureOutpath <- paste0(outputDir, name, ".png")
      cat("Figure output filepath:", figureOutpath, "\n")
      png(filename = figureOutpath, width = plotWidth, height = plotHeight, res = plotDPI)
      
      #### Generate the track plot ####
      plotTracks(list(track.ideogram, track.axis,
                      track.bigwig.prog_A, track.bigwig.prog_B,
                      track.bigwig.luf_A, track.bigwig.luf_B,
                      track.bigwig.baz2b_A, track.bigwig.baz2b_B,
                      track.UCSCgenes.symbols),
                 from = plot.from, to = plot.to, sizes = plot.sizes,
                 main = plot.title, cex.main = 1.25, transcriptAnnotation = "symbol")
      
      ##
      dev.off()
      
    }, finally = {
      next
    })
    
  } # end for (a in 1:numItems)
} # end generatePlotsBigwig function

#### Generate BAM plots ####
generateBamPlots <- function(chrList, startList, endList, nameList, outputDir, plotWidth, plotHeight, plotDPI){
  
  ## Report
  numItems <- length(chrList)
  cat("Generating BAM plots", "\n")
  cat("Number of items:", numItems, "\n")
  cat("Output directory for figures:", outputDir, "\n")
  cat("Plot width:", plotWidth, "\n")
  cat("Plot height:", plotHeight, "\n")
  cat("Plot DPI:", plotDPI, "\n")
  
  ## Loop
  for (a in 1:numItems){
    
    tryCatch({
      
      ## Get current data
      chr <- chrList[a]
      start <- startList[a]
      end <- endList[a]
      name <- nameList[a]
      
      ##
      cat("Generating plot number:", a, "\n")
      cat("Current chromosome:", chr, "\n")
      cat("Current start:", start, "\n")
      cat("Current end:", end, "\n")
      cat("Current gene:", name, "\n")
      
      ## Set paramters for track initiation
      cat("Setting track parameters", "\n")
      setParameters("hg38", chr, start, end, "*",
                    TRUE, 0.5, c(0,50), "coverage",
                    "histogram", c(0,20), 1, "black")
      
      ##
      cat("Initializing tracks", "\n")
      track.ideogram <- IdeogramTrack(genome = genome, chromosome = chromosome, start = start, end = end)
      track.axis <- GenomeAxisTrack(genome = genome, chromosome = chromosome, start = start, end = end)
      
      ## The gene annotations track with gene symbols
      ucscGenes <- UcscTrack(genome=genome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                             chromosome=chromosome, rstarts = "exonStarts", rends = "exonEnds",
                             gene = "name", symbol = 'name', transcript = "name",
                             strand = "strand", stacking = 'pack', showID = T, geneSymbol = T,
                             showTitle = T, name = "GENE", fontcolor.title = "white",
                             lineheight=1, rotation.group=90, just.group="left", fontsize.group=12)
      z <- ranges(ucscGenes)
      mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
      ucscGenes2 <- ucscGenes
      ranges(ucscGenes2) <- z
      track.UCSCgenes.symbols <- ucscGenes2
      
      #### Initiate the BAM tracks ####
      cat("Initializing BAM tracks", "\n")
      track.bam.prog_A <- makeBamTrack(path.bam.progenitor_A, "PROG_A", "black", "black", "blue")
      track.bam.prog_B <- makeBamTrack(path.bam.progenitor_B, "PRoG_B", "black", "black", "blue")
      track.bam.luf_A <- makeBamTrack(path.bam.luciferase_A, "LUF_A", "black", "black", "yellow")
      track.bam.luf_B <- makeBamTrack(path.bam.luciferase_B, "LUF_B", "black", "black", "yellow")
      track.bam.baz2b_A <- makeBamTrack(path.bam.baz2b_A, "BAZ2B_A", "black", "black", "red")
      track.bam.baz2b_B <- makeBamTrack(path.bam.baz2b_B, "BAZ2B_B", "black", "black", "red")
      
      ## Set plotting parameters
      cat("Setting plotting parameters", "\n")
      setPlotParameters(c(0.5,0.5,1,1,1,1,1,1,2), 5000, 5000, name)
      
      ## Generate the plot
      cat("Generating plot", "\n")
      figureOutpath <- paste0(outputDir, name, ".png")
      cat("Figure output filepath:", figureOutpath, "\n")
      png(filename = figureOutpath, width = plotWidth, height = plotHeight, res = plotDPI)
      
      #### Generate the track plot ####
      plotTracks(list(track.ideogram, track.axis,
                      track.bam.prog_A, track.bam.prog_B,
                      track.bam.luf_A, track.bam.luf_B,
                      track.bam.baz2b_A, track.bam.baz2b_B,
                      track.UCSCgenes.symbols),
                 from = plot.from, to = plot.to, sizes = plot.sizes,
                 main = plot.title, cex.main = 1.25, transcriptAnnotation = "symbol")
      
      ##
      dev.off()
      
    }, finally = {
      next
    })
  } # end for (a in 1:numItems)
} # end generatePlots function

#### Generate BAM tracks ####
makeBamTrack <- function(bamPath, trackName, titleColor, covOutlineColor, covBodyColor){
  
  ####
  tempTrack <- AlignmentsTrack(range = bamPath,
                               genome = genome, chromosome = chromosome,
                               start = start, end = end,
                               strand = strand,
                               cex = at.cex,
                               type = at.type,
                               showTitle = at.showTitle,
                               ylim = at.ylim,
                               name = trackName, # title of the track
                               fontcolor.title = titleColor, # color of the title
                               col.coverage = covOutlineColor, # the outline color of the coverage lines
                               fill.coverage = covBodyColor) # the fill color of the coverage body
  
  ####
  return(tempTrack)
  
} # end makeBamTrack function













