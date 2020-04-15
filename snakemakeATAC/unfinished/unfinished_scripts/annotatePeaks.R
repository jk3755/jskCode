## Set snakemake variables
cat("Setting snakemake variables", "\n")
bedFile <- snakemake@input[[1]]
functionSourcePath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]

## Source functions ##
cat("Loading functions from path:", functionSourcePath, "\n")
source(functionSourcePath)

##
filePathBase <- gsub("annotate.peaks.complete", "", outPath)

## Profile of peaks in TSS regions
cat("Generating TSS peak profile", "\n")
TSSprofilePath <- paste0(filePathBase, ".TSSprofile.svg")
plotTSSpeakProfile(bedFile, TSSprofilePath)

## Average profile of ChIP peaks binding to TSS region
cat("Generating average TSS peak profile", "\n")
avgTSSprofilePath <- paste0(filePathBase, ".avgTSSprofile.svg")
plotAvgTSSpeakProfile(bedFile, avgTSSprofilePath)

## Peak annotations
cat("Generating peak annotations", "\n")
annoPlot1 <- paste0(filePathBase, ".annoplot1.svg")
annoPlot2 <- paste0(filePathBase, ".annoplot2.svg")
annoPlot3 <- paste0(filePathBase, ".annoplot3.svg")
annoPlot4 <- paste0(filePathBase, ".annoplot4.svg")
generatePeakAnnotationPlots(bedFile, annoPlot1, annoPlot2, annoPlot3, annoPlot4)

## Finish up
cat("Touching snakemake flag file", "\n")
file.create(outPath)