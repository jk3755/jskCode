
#### Load the counts (tpm) ####
ccleCountsPath <- "C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\LNCaP\\counts\\CCLE_RNAseq_rsem_genes_tpm_20180929.txt"
ccleCountsTPM <- read.table(ccleCountsPath, sep="\t", header = TRUE)

##
ensembl <- as.character(ccleCountsTPM[,1])
symb <- getGeneSymbolsFromEnsemblId(ensembl)
paadCountsSymbol <- paadCountsEntrez
rownames(paadCountsSymbol) <- symb

msviperAndPlot <- function( signature , interactomeObj , path = "reports" , filetag = "marina-plot" ,
                            mrs = 25 , returnMarinaObj = FALSE , nullModel = NULL , purgeSignature = FALSE )
{
  if( purgeSignature )
  {
    signature <- purgeSignature( signature , interactomeObj )
    # nullModel.purged <- subset( nullModel , rownames(nullModel) %in% names(signature.purged) )
    message(" >>> Signature purged !!!")
    
  }
  # protein_activity <- msviper( ges = signature.purged , regulon = interactomeObj , verbose = T , cores = parallel::detectCores() )
  if (!is.null(nullModel))
  {
    message("*** Usign a Null Model ***")
    protein_activity <- msviper( signature , regulon = interactomeObj , nullmodel = nullModel )
  }
  else
  {
    protein_activity <- msviper( signature , regulon = interactomeObj )
  }
  
  filename <- file.path( path , paste0( "protein_activity-" , filetag , ".pdf" ) )
  pdf( file = filename )
  plot( protein_activity , mrs = mrs )
  dev.off()
  
  if (returnMarinaObj)
    return(protein_activity)
}

msviperAndPlot(paadCounts, regulonpradSYMBOLS, path = getwd())

##########################################################################################################################
ccleCountsPath <- "C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\LNCaP\\counts\\cclecounts1.rda"
load(ccleCountsPath)

# matrix is called ccle
colnames(ccle)

sampleNames <- c("LNCAPCLONEFGC_PROSTATE", "22RV1_PROSTATE", "DU145_PROSTATE", "NCIH660_PROSTATE", "PC3_PROSTATE", "MDAPCA2B_PROSTATE", "VCAP_PROSTATE")

sampleNames <- c("LNCAPCLONEFGC_PROSTATE")
names <- colnames(ccle)

## indices
lncap <- 739
rv1 <- 734
du145 <- 735
h660 <- 736

##
idx <- c(739,734,735,736)

lncapMatrix <- ccle[,idx]
refMatrix <- ccle[,-idx]
colnames(lncapMatrix) <- "LNCAP"

# seems like need to run viper signature with more than one sample
x <- viperSignature(lncapMatrix, refMatrix, method="zscore")

y <- viper(x, regulonprad, method="none")



Get RPKM normalized counts matrix from CCLE
```{r}

## Read in LNCAP WT expression counts (RPKM) from GCT file
countsGCTpath <- "C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\LNCaP\\counts\\CCLE_RNAseq_genes_rpkm_20180929.gct"
gctCCLE <- read.gct(countsGCTpath)
countsCCLE <- gctCCLE@assayData[["exprs"]]

##
sampleName1 <- "LNCAPCLONEFGC_PROSTATE"
sampleName2 <- "22RV1_PROSTATE"
sampleName3 <- "DU145_PROSTATE"
sampleName4 <- "NCIH660_PROSTATE"
sampleName5 <- "PC3_PROSTATE"
sampleName6 <- "MDAPCA2B_PROSTATE"
sampleName7 <- "VCAP_PROSTATE"

sampleNames <- c("LNCAPCLONEFGC_PROSTATE", "22RV1_PROSTATE", "DU145_PROSTATE", "NCIH660_PROSTATE", "PC3_PROSTATE", "MDAPCA2B_PROSTATE", "VCAP_PROSTATE")

##
eset1 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName1]
eset2 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName2]
eset3 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName3]
eset4 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName4]
eset5 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName5]
eset6 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName6]
eset7 = gctCCLE[, sampleNames(gctCCLE) %in% sampleName7]

##
counts1 <- eset1@assayData[["exprs"]]
counts2 <- eset2@assayData[["exprs"]]
counts3 <- eset3@assayData[["exprs"]]
counts4 <- eset4@assayData[["exprs"]]
counts5 <- eset5@assayData[["exprs"]]
counts6 <- eset6@assayData[["exprs"]]
counts7 <- eset7@assayData[["exprs"]]

## Change to gene symbols
geneSymbol <- gctCCLE@featureData@data[["Description"]]

##
prostateCountMatrix.RPKM <- matrix(data = NA, ncol = 7, nrow = 56202)
rownames(prostateCountMatrix.RPKM) <- geneSymbol
colnames(prostateCountMatrix.RPKM) <- sampleNames

##
prostateCountMatrix.RPKM[,1] <- counts1
prostateCountMatrix.RPKM[,2] <- counts2
prostateCountMatrix.RPKM[,3] <- counts3
prostateCountMatrix.RPKM[,4] <- counts4
prostateCountMatrix.RPKM[,5] <- counts5
prostateCountMatrix.RPKM[,6] <- counts6
prostateCountMatrix.RPKM[,7] <- counts7

```

Process the RPKM matrix (simpleviper)
```{r}

# 1.1 Get a list of all targets of every regulator in the network.
tmp <- c(names(regulonpaad), unlist(lapply(regulonpaad, function(x) names(x$tfmode)), use.names = FALSE))

# 1.2 Remove expression data of all the genes which arent in the list identified above.
prostateCountMatrix.RPKM <- prostateCountMatrix.RPKM[rownames(prostateCountMatrix.RPKM) %in% unique(tmp),]

# 1.3 Scale the expression matrix on rows [i.e. row = (row - mean(row)) / sd(row) ].
tt <- t(scale(t(prostateCountMatrix.RPKM)))

# 2.1 Remove targets in the regulon which are not in the rownames of the expression matrix.
regulonpaad <- lapply(regulonpaad, function(x, genes) {
  filtro <- names(x$tfmode) %in% genes
  x$tfmode <- x$tfmode[filtro]
  if (length(x$likelihood) == length(filtro))
    x$likelihood <- x$likelihood[filtro]
  return(x)
}, genes = rownames(prostateCountMatrix.RPKM))

# 2.2 Define minimum regulon size for filtering (default is 20).
minsize <- 20

# 2.3 Remove regulators with a regulon size below the minsize parameter.
regulonpaad <- regulonpaad[sapply(regulonpaad, function(x) length(x$tfmode)) >= minsize]

# 1.1 Get a list of all targets of every regulator in the network.
targets <- unique(unlist(lapply(regulonpaad, function(x) names(x$tfmode)), use.names = FALSE))

# 1.2 Create the Mode of Regulation matrix from the regulon object.
mor <- sapply(regulonpaad, function(x, genes) {
  return(x$tfmode[match(genes, names(x$tfmode))])
}, genes = targets)

# 1.2 Create the Weights matrix from the regulon object.
wts <- sapply(regulonpaad, function(x, genes) {
  tmp <- x$likelihood[match(genes, names(x$tfmode))]
  tmp[is.na(match(genes, names(x$tfmode)))] <- NA
  return(tmp/max(tmp, na.rm = T))
}, genes = targets)

# 1.3 For each regulator, assign values of 0 to genes which are not listed as its targets.
mor[is.na(mor)] <- 0
wts[is.na(wts)] <- 0

# 1.4 Scale the columns of the Weights matrix to the sum of the weights.
wtss <- scale(wts, center = FALSE, scale = colSums(wts))

# 2.1 Calculate the T2 rank matrix from the expression dataset.
t2 <- apply(tt, 2, rank)/(nrow(tt) + 1)

# This line of code is necessary to match the order of genes
# for the subsequent matrix multiplication steps.
pos <- match(targets, rownames(tt))

# 2.2 Transform T2 ranks to Gaussian values.
t2q <- qnorm(filterRowMatrix(t2, pos))

# 2.3 Matrix multiplication.
sum1 <- t(mor * wtss) %*% t2q

# 3.1 Calculate the T1 Rank matrix from the T2 Rank matrix.
t1 <- abs(t2 - 0.5) * 2
t1 <- t1 + (1 - max(t1))/2

# 3.2 Transform T1 ranks to Gaussian values.
t1q <- qnorm(filterRowMatrix(t1, pos))

# 3.3 Matrix multiplication.
sum2 <- t((1 - abs(mor)) * wtss) %*% t1q

# 4.1 Extract the signs of the Two-tail enrichment scores
ss <- sign(sum1)
ss[ss == 0] <- 1

# 4.2 Combine the Two-tail and One-tail enrichment score matrices.
sum3 <- (abs(sum1) + sum2 * (sum2 > 0)) * ss

# 5.1 For each regulator, calculate an index proportional to the likelihood value
# of all its regulatory interactions.
lwt <- sqrt(colSums(wts^2))

# 5.2 Adjust 3T enrichment scores proportionally to the weights calculated above.
nes <- sum3 * lwt

```

.libPaths("C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\viper\\lib")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("mixtools")
BiocManager::install("bcellViper")
BiocManager::install("viper")
BiocManager::install("aracne.networks")
library(viper)
library(aracne.networks)


data(regulonprad)
regulonPRADentrez <- regulonprad
regulonPRADsymbol <- entrez2gene.regulon(regulonPRADentrez)

ccleCountsPath <- "C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\LNCaP\\counts\\cclecounts.rda"
load(ccleCountsPath)
paadCountsEntrez <- cclecounts[["prostate_bat1"]]
entrez <- rownames(paadCountsEntrez)
symb <- entrez2geneSymbol.old(entrez)
paadCountsSymbol <- paadCountsEntrez
rownames(paadCountsSymbol) <- symb


ccleCountsPath <- "C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\LNCaP\\counts\\cclecounts1.rda"
load(ccleCountsPath)

# matrix is called ccle
colnames(ccle)

sampleNames <- c("LNCAPCLONEFGC_PROSTATE", "22RV1_PROSTATE", "DU145_PROSTATE", "NCIH660_PROSTATE", "PC3_PROSTATE", "MDAPCA2B_PROSTATE", "VCAP_PROSTATE")
names <- colnames(ccle)

## indices
lncap <- 739
rv1 <- 734
du145 <- 735
h660 <- 736

##
idx <- c(739,734,735,736)

lncapMatrix <- ccle[,idx]
refMatrix <- ccle[,-idx]
#colnames(lncapMatrix) <- "LNCAP"

# seems like need to run viper signature with more than one sample
x <- viperSignature(lncapMatrix, refMatrix, method="zscore")
x2 <- x[["signature"]]

y <- viper(x, regulon = regulonPRADentrez, method="none", verbose = T)


##
nullmodel <- ttestNull(lncapMatrix, refMatrix, per = 1000, repos = TRUE, verbose = TRUE)


####
data(regulonprad)
regulonPRADentrez <- regulonprad
ccleCountsPath <- "C:\\Users\\jsk33\\OneDrive\\R\\ATACseqAnalysis\\LNCaP\\counts\\cclecounts.rda"
load(ccleCountsPath)
idx <- c(739,734,735,736)
lncapMatrix <- ccle[,idx]
refMatrix <- ccle[,-idx]
x <- viperSignature(lncapMatrix, refMatrix, method="zscore")
x$nullmodel[which(is.na(x$nullmodel))] <- 0
x$nullmodel[which(is.infinite(x$nullmodel))] <- 0
x$signature[which(is.infinite(x$signature))] <- 0
y <- viper(x, regulonPRADentrez, method="none")


idx <- which(is.na(x[["signature"]]) == TRUE)
x[["signature"]][idx] <- 0

idx <- which(is.na(x[["nullmodel"]]) == TRUE)
x[["signature"]][idx] <- 0

idx <- which(is.na(x[["signature"]]) == TRUE)
x[["signature"]][idx] <- 0


x$signature[which(is.infinite(x$signature))] <- 0
