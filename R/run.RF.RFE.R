##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Running Recursive feature elimination (RFE) with a Random forest (RF) model on the training dataset
## Input: Training data frame (trainData), number of features to be tested (subsets)
## Output: Trained classification model
## Date: 28-9-2018
##------------------------------------------

run.RF.RFE <- function(trainData, subsets = c(2:30)){
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(caret)

  #### Train the RFE-RF classification model on the training set using 5 x 10-fold cross validation
  cat("Running Recursive feature elimination (RFE) with a Random forest (RF) model using 5 times 10-fold CV..\n")
  ctrl <- rfeControl(functions = rfFuncs,
                     method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     verbose = TRUE,
                     saveDetails = TRUE,
                     returnResamp = "all")

  trainctrl <- trainControl(classProbs= TRUE,
                            verboseIter = FALSE,
                            returnResamp = "final",
                            returnData = TRUE)
  set.seed(12)
  rfe <- rfe(class ~.,
             data = trainData,
             sizes = subsets,
             method="rf",
             rfeControl = ctrl,
             trControl = trainctrl,
             preProc = c("center", "scale"))
  
  return(rfe)
}