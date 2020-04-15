
#### Plot the correlation between bam files ####


####
makeGenomicWindows <- function(binsize){
  bins = tileGenome(seqinfo(Hsapiens),
                    tilewidth = binsize,
                    cut.last.tile.in.chrom = TRUE)
  bins <- cleanGRanges(bins)
  return(bins)
}

#### Remove PhiX Reads from a set of paired fastq files ####
pairedRemovePhiX <- function(inputFilePaths, outputFilePaths){
  
  ##
  cat("Running pairedRemovePhiX function", "\n")
  
  ## Library
  cat("Loading libraries", "\n")
  #BiocManager::install("dada2")
  suppressMessages(library(dada2))
  
  ##
  cat("Read 1 input:", inputFilePaths[1], "\n")
  cat("Read 2 input:", inputFilePaths[2], "\n")
  cat("Read 1 output:", outputFilePaths[1], "\n")
  cat("Read 2 output:", outputFilePaths[2], "\n")
  
  ## Remove the phiX and save the output files
  cat("Removing PhiX reads and saving output files", "\n")
  fastqClean <- fastqPairedFilter(inputFilePaths, outputFilePaths, rm.phix = TRUE)
  
} # end pairedRemovePhiX

####
getTotalReadsFromBam <- function(inputBam){
  
  ## Load required libraries
  suppressMessages(library(Rsamtools))
  
  ## Count total reads in library and calculate average signal per bp
  totalReads <- countBam(inputBam)
  totalReads <- as.numeric(totalReads[6])
  
  ## Return
  return(totalReads)
  
} # end getTotalReadsFromBam

####
removeOutlierWidths <- function(inputWidths){
  outlierWidths <- boxplot(inputWidths, plot = FALSE)$out
  outlierIdx <- which(inputWidths %in% outlierWidths)
  inputWidths <- inputWidths[-outlierIdx]
  return(inputWidths)}

#### Retrieve underlying base sequence from a GRanges object
getUnderlyingSequence <- function(grInput){
  
  ## Report
  cat("Running getUnderlyingSequence function", "\n")
  
  ## Load required libraries
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  genome <- Hsapiens
  
  ## Get the underlying sequence
  siteSequences <- Views(genome, grInput)
  siteSequences <- as.character(siteSequences)
  
  ## Return the character sequences
  return(siteSequences)}

#### Create a fasta sequence file of hg38
generateFastaHg38 <- function(outputPath){
  
  ##
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  genome <- Hsapiens
  
  ## 
  export(genome, outputPath, format = "fasta")}

#### Return a DF with hg38 chromosome and their lengths
getChrLengthsHg38 <- function(){
  
  chr1len <- 248956422
  chr2len <- 242193529
  chr3len <- 198295559
  chr4len <- 190214555
  chr5len <- 	181538259
  chr6len <- 170805979
  chr7len <- 159345973
  chr8len <- 145138636
  chr9len <- 	138394717
  chr10len <- 133797422
  chr11len <- 135086622
  chr12len <- 133275309
  chr13len <- 114364328
  chr14len <- 107043718
  chr15len <- 101991189
  chr16len <- 90338345
  chr17len <- 83257441
  chr18len <- 80373285
  chr19len <- 58617616
  chr20len <- 64444167
  chr21len <- 46709983
  chr22len <- 50818468
  chrXlen <- 156040895
  chrYlen <- 57227415
  #
  lengths <- c(chr1len, chr2len, chr3len, chr4len, chr5len, chr6len, chr7len, chr8len, chr9len, chr10len,
               chr11len, chr12len, chr13len, chr14len, chr15len, chr16len, chr17len, chr18len, chr19len, chr20len,
               chr21len, chr22len, chrXlen, chrYlen)
  #
  id = factor(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),
              levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                         "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  #
  df <- data.frame(chr = id, length = lengths, stringsAsFactors = FALSE)
  #
  return(df)
}

#### Calculate the total number of reads based on a peak file
getTotalReadsInPeaks <- function(inputBam, peaksPath){
  
  ## Load required libraries
  suppressMessages(library(Rsamtools))
  
  ## Import the peaks
  accessibilityPeaks <- importMACS2narrowPeak(peaksPath)
  accessibilityPeaks <- cleanGRanges(accessibilityPeaks)
  
  ##
  param = ScanBamParam(which = accessibilityPeaks)
  totalReadsInPeaks <- countBam(inputBam, param = param)
  totalReadsInPeaks <- as.numeric(sum(totalReadsInPeaks[,6]))
  
  ##
  return(totalReadsInPeaks)
  
} # end getTotalReadsInPeaks

#### Calculate total reads in a bam file in a range provided by Granges object
getTotalReadsInGRanges <- function(inputBam, inputGR){
  
  ## Load required libraries
  suppressMessages(library(Rsamtools))
  
  ##
  param = ScanBamParam(which = inputGR)
  
  ##
  totalReadsInRange <- countBam(inputBam, param = param)
  
  ##
  totalReadsInRange <- as.numeric(sum(totalReadsInRange[,6]))
  
  ##
  return(totalReadsInRange)
  
} # end getTotalReadsInGRanges

#### Convert genrich peaks to HOMER input peaks ####
convertGenrichPeakToHomerFormat <- function(genrichPath, homerPath, sampleName){
  
  ## read in the genrich file
  genrichPeaks <- read.table(genrichPath)
  
  ## add the column names
  ## ucsc format https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  cnames <- c("chr", "start", "end", "name", "score", "strand", "signal", "pvalue", "qvalue", "summit")
  colnames(genrichPeaks) <- cnames
  
  ## assign unique sample name to peak names
  genrichPeaks$name <- paste0(sampleName, "_", genrichPeaks$name)
  
  ## assign strand
  genrichPeaks$strand <- 0
  
  ## retrieve required data for homer peaks
  homerPeaks <- data.frame(name = genrichPeaks$name,
                           chr = genrichPeaks$chr,
                           start = genrichPeaks$start,
                           end = genrichPeaks$end,
                           strand = genrichPeaks$strand)
  
  ## write to file
  write.table(homerPeaks,
              file = homerPath,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  
} # end convertGenrichPeakToHomerFormat function

#### Convert genrich peaks to HOMER input peaks ####
writeGRangesToHomerFormat <- function(ranges, outputPath){
  
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
              file = outputPath,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
}












