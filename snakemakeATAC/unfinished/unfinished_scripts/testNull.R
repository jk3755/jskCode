
##
BiocManager::install("seqbias")
library(seqbias)

##
bindingSites <- getAllBindingSites("MNX1", pwmScanScore = "95%")
bindingSite <- bindingSites[1]

##
seqbiasPath <- "C:\\Users\\Jordan\\OneDrive\\labData\\seqbias\\human\\SRR1554094.hg38.chr1.yaml"
fastaPath <- "C:\\Users\\Jordan\\OneDrive\\labData\\seqbias\\human\\chr1.fa"
totalSignal <- 500
upstream <- 100
downstream <- 100
iterations <- 1000

##
testCorr <- generateNullFootprintSeqbiasCorrected(totalSignal, bindingSite, seqbiasPath, fastaPath)

##
test <- generateNullFootprint(totalSignal, bindingSite)
