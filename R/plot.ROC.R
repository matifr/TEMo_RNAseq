##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Draw ROC curves on the training and testing datasets
## Input: Trained model (trainModel), predictions on the testData set (testModel), class of the testData set (testDataClass)
## Output: ROC curves plot of the model on the training and testing datasets
## Date: 5-10-2018
##------------------------------------------

plot.ROC = function(trainModel, testModel, testDataClass, Title = "ROC Curve"){
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  library(ROCR)

  pred = prediction(predictions = trainModel$pred$Normal, labels = trainModel$pred$obs, label.ordering = rev(levels(c("Cancer", "Normal"))))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, main=Title, col ="#1B9E77",lwd=2, cex.axis=10)
  abline(0, 1, col="black", lty=2)
  
  auc <- performance(pred,"auc")
  auc <- unlist(slot(auc, "y.values"))
  
  auc<-round(auc, digits = 2) * 100
  auct <- paste(c("Training (AUC: "),auc,c("%)"))
  
  pred <- prediction(predictions = testModel$Normal, labels = testDataClass, label.ordering = rev(levels(c("Cancer", "Normal"))))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, col ="#7570B3",lwd=2, add=TRUE)
  
  auc1 <- performance(pred,"auc")
  auc1 <- unlist(slot(auc1, "y.values"))
  
  auc1<-round(auc1, digits = 2) * 100
  auct1 <- paste(c("Testing (AUC: "),auc1,c("%)"))
  
  legend(0.5,0.3, legend=c(c(auct, auct1)),cex=1,text.col=c('#1B9E77', "#7570B3"), col=c("#1B9E77", "#7570B3") ,bty = "n", text.font = 3)
}