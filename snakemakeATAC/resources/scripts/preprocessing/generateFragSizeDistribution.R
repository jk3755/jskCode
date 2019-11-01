#### Set snakemake variables ####
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
functionSourcePath <- snakemake@input[[3]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
refGenome <- snakemake@wildcards[["refgenome"]]

## Source functions ##
cat("Loading functions from path:", functionSourcePath, "\n")
source(functionSourcePath)

#### Output paths ####
svgOut1 <- snakemake@output[[1]]
svgOut2 <- snakemake@output[[2]]
svgOut3 <- snakemake@output[[3]]
svgOut4 <- snakemake@output[[4]]

#### Report ####
cat("Bam filepath:", bamPath, "\n")
cat("Bai filepath:", baiPath, "\n")
cat("Sample name:", sampleName, "\n")
cat("Sample rep:", sampleRep, "\n")
cat("Reference genome:", refGenome, "\n")
cat("Output for svg plot 1:", svgOut1, "\n")
cat("Output for svg plot 2:", svgOut2, "\n")
cat("Output for svg plot 3:", svgOut3, "\n")
cat("Output for svg plot 4:", svgOut4, "\n")

#### Make the plots ####
plotFragmentSizeDistribution(bamPath, svgOut1, svgOut2, svgOut3, svgOut4)
