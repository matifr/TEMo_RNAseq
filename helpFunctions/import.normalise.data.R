##------------------------------------------
## Author: Matina Fragkogianni
## Date: 28-9-2018
## This function imports the gene expression data and create a normalized matrix of gene expression for downstream analysis
## 
## @return A matrix containing the normalised expression values
## 
##------------------------------------------

Import.normalise.data <- function(){
 
  ##------------------------------------------
  # Source files
  ##------------------------------------------
   source(file = "helpFunctions/remove.dots.R")
  
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(limma)
  require(edgeR)
  require(ggplot2)
  require(ggrepel)
  
  #### Read sample information and raw counts
  pData <- read.delim(file = "phenoData/ALL_pData.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  pData <- pData[which(pData$Batch_info %in% c("4","5","6")),]
  pData <- pData[which(pData$Condition %in% c("Breast_cancer","Normal")),]
  pData$Condition2 <- ifelse(pData$Condition == "Normal","Normal","Cancer")
  
  cat("Reading files in..")
  d <- readDGE(pData$Filenames, path = "Counts/", columns = c(1,2), pData$Condition2)
  d$samples$batch = pData$Batch_info
  
  #### Filter lowly expressed genes and recalculate library size
  cat("Filtering lowly expressed genes..")
  keep <- rowSums(cpm(d) > 1) >= 32
  d <- d[keep,]
  print(table(keep))
  colnames(d) <- pData$SampleName
  d$samples$libSize <- colSums(d$counts)
  
  #### Upperquartile normalization and log2 transformation
  cat("Upperquartile normalization and log2 transformation..")
  d <- calcNormFactors(d, method = "upperquartile")
  lcpm <- cpm(d,log = TRUE,prior.count = 1)
  
  #### Remove unecessery dots from the ensemble IDs
  rownames(lcpm) <- remove.dots(lcpm)
  
  return(list("lcpm" = lcpm,"dge" = d))
}
