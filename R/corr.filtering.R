##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Dimensionality reduction based on feature correlation on the training data
## Input: Training data matrix (trainData)
## Output: List containing the reduced gene expression matrix and dataFrame containing sample class information
## Date: 1-10-2018
##------------------------------------------

corr.filtering = function(trainData){
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(feseR)
  require(dplyr)
  
  
  tData = trainData %>% dplyr::select(-c("class"))
  cT =  trainData %>% select(c("class"))
  cT$class= ifelse(cT$class == "Normal","0","1")
  tData = as.matrix(tData)
  cT$class = as.numeric(cT$class)
  cT = as.matrix(cT$class)
  reduced.tData = filter.corr(features = tData, class = cT, mincorr = 0.4)
  dim(reduced.tData)
  
  
  cT = as.data.frame(cT)
  names(cT) = "class"
  cT$class = ifelse(cT$class == "0","Normal","Cancer")
  cT$class = as.factor(cT$class)
  
  return(list("reduced.tData" = reduced.tData, "class" = cT))
  
}