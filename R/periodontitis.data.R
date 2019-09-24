##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Import and processing of the public dataset of isolated monocytes from periodontitis samples
## Input: Gene signature selected during model training (predictors)
## Output: Dataframe containing gene expression and sample class information
## Date: 1-10-2018
##------------------------------------------
  
periodontitis.data <- function(predictors){
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
  nrmData = cpm(d,log=TRUE,prior.count=1)
  
  
  rownames(nrmData)=remove.dots(nrmData)
  cat(paste("Number of predictors in the data: ",sum(rownames(nrmData) %in% predictors)))
  
  #### Create data frame for the data
  periodTestData = t(nrmData)
  periodTestDataClass = periodSampleConditions
  periodTestData = data.frame(periodTestData, class = periodTestDataClass)
  
  return(periodTestData)
}