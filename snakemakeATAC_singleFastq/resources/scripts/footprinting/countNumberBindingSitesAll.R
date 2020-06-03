
####
suppressMessages(library(stringr))

#### Generate the input directory path ####
inputDirectory <- "C:\\Users\\Jordan\\OneDrive\\LNCAP_count_binding_sites\\data"

#### Get the input file paths ####
inputFileList <- list.files(path = inputDirectory, full.names = TRUE)
numInputFiles <- length(inputFileList)

#### 
geneNames <- c()
numBindingSites <- c()

####
for (a in 1:numInputFiles){
  
  ##
  tempPath <- inputFileList[a]
  
  ## Get gene name
  tempGene <- str_replace(tempPath, "C:\\Users\\Jordan\\OneDrive\\LNCAP_count_binding_sites\\data", "")
  tempGene <- str_replace(tempGene, ".RData", "")
  cat(tempGene, "\n")
  
  
  ##
  load(inputFileList[a])
  
  ##
  
  
}

####
sitesInfo <- data.frame(gene = geneNames,
                        numBindingSites = numBindingSites)

####
outPath <- "C:\\Users\\Jordan\\OneDrive\\LNCAP_count_binding_sites\\sitesInfo.RDS"
saveRDS(sitesInfo, file = dataOutPath)

