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
  
  #### Load the DGE object containing the periodontitis dataset
  load(file = "Data/periodontitis_dge_obj.Rdata")
  
  d <- calcNormFactors(d, method = "upperquartile") 
  nrm_data = cpm(d,log=TRUE,prior.count=1)
  
  
  rownames(nrm_data)=remove.dots(nrm_data)
  cat(paste("Number of predictors in the data: ",sum(rownames(nrm_data) %in% predictors)))
  
  period_testData = t(nrm_data)
  period_testData_class = period_sampleConditions
  period_testData = data.frame(period_testData, class = period_testData_class)
  
  return(period_testData)
}