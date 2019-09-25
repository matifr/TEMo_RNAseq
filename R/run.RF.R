##------------------------------------------
## Author: Matina Fragkogianni
## Date: 5-11-2018
## 
## A function to train a random forest model
## @param trainData a dataframe , 
## @return The trained random forest model
## 
##------------------------------------------

run.RF <- function(trainData){
  
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