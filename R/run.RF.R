##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Trains  random forest model
## Input: Training dataset (TrainData), 
## Output: The trained random forest model
## Date: 5-11-2018
##------------------------------------------

run.RF = function(trainData){
  
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(caret)

  ctrl <- trainControl(method = "repeatedcv", number = 10,
                       repeats = 5,
                       classProbs = TRUE,
                       verboseIter = FALSE, 
                       savePredictions = TRUE,
                       returnResamp = "final")

  # Train the model using random forest
  set.seed(12)
  model_rf <- train(class ~. ,
                    trainData,
                    method = "rf",
                    #metric = "ROC",
                    ntree = 500,
                    trControl = ctrl, 
                    importance = TRUE,
                    preProc=c("center", "scale"),
                    tuneGrid = expand.grid(.mtry=c(4))
                    )

  return(model_rf)
}