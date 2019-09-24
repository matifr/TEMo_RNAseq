##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Import and processing of the public dataset of isolated monocytes from periodontitis samples
## Input: Gene signature selected during model training (predictors)
## Output: Dataframe containing gene expression and sample class information
## Date: 1-10-2018
##------------------------------------------
  
periodontitis.data = function(predictors){
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(edgeR)
  require(caret)
  
  ##------------------------------------------
  # Source files
  ##------------------------------------------
  source(file = "Functions/remove.dots.R")
  
  sampleFiles = list.files(path = "Input/public_periodontitis/", pattern = "*.count")
  period_sampleConditions = c(rep("Cancer",5), rep("Normal",5))
  d = readDGE(sampleFiles, path = "Input/public_periodontitis/", columns=c(1,2), period_sampleConditions)
  
  d <- calcNormFactors(d, method = "upperquartile") 
  nrm_data = cpm(d,log=TRUE,prior.count=1)
  
  
  rownames(nrm_data)=remove.dots(nrm_data)
  cat(paste("Number of predictors in the data: ",sum(rownames(nrm_data) %in% predictors)))
  
  period_testData = t(nrm_data)
  period_testData_class = period_sampleConditions
  period_testData = data.frame(period_testData, class = period_testData_class)
  
  return(period_testData)
}