##------------------------------------------
## Author: Matina Fragkogianni
## Date: 7-11-2018
##
## This is the main script that performs model training and classification using monocyte gene expression
## 
##------------------------------------------

##------------------------------------------
# Load libraries
##------------------------------------------
library(pROC)
library(caret)
library(dplyr)
library(ROCR)

##------------------------------------------
# Source files
##------------------------------------------
source(file = "helpFunctions/import.normalise.data.R")
source(file = "helpFunctions/run.ComBat.R")
source(file = "helpFunctions/my.pca.R")
source(file = "helpFunctions/corr.filtering.R")
source(file = "helpFunctions/periodontitis.data.R")
source(file = "helpFunctions/ensembl.to.GS.R")
source(file = "helpFunctions/plot.ROC.R")
source(file = "helpFunctions/run.RF.RFE.R")

##------------------------------------------
# Data import and model training and testing  
##------------------------------------------

#### Import and normalise gene expression monocyte data
dlist <- Import.normalise.data()

#### PCA plot before batch effect correction
my.pca(dlist$lcpm, dge = dlist$dge)

#### Correct data for batch effects using Combat 
combatEdata <- run.ComBat(dge = dlist$dge, mat = dlist$lcpm)

#### PCA plot after batch effect correction
my.pca(combatEdata, dge = dlist$dge)

data <- combatEdata
class <- as.factor(dlist$dge$samples$group)
class <- relevel(class, "Cancer")
data <- t(data)
dim(data)
data_fr <- data.frame(data, class = class)

#### Split dataset into training (70%) and testing (30%)
set.seed(12)

# Step 1: Get row numbers for the training data
trainRowNumbers <- createDataPartition(data_fr$class, p = 0.7, list=FALSE)

# Step 2: Create the training  dataset
trainData <- data_fr[trainRowNumbers,]
dim(trainData)

# Step 3: Create the test dataset
testData <- data_fr[-trainRowNumbers,]
dim(testData)


#### Filtering based on linear correlation
reducedDataList <- corr.filtering(trainData = trainData)
reducedtData <- reducedDataList$reducedtData
dim(reducedtData)
reducedtData.class <- reducedDataList$class
reducedtData <- as.data.frame(reducedtData)
reducedtData$class <- reducedtData.class$class

#### Classification and feature selection using RFE-RF on the training set
rfCorrRFE <- run.RF.RFE(trainData = trainData)
confusionMatrix(rfCorrRFE,positive = "Cancer")
ensembl.to.GS(rfCorrRFE$optVariables)
plot(rfCorrRFE, type = c("g", "o"), xlim = c(0:31), ylim = c(0.7,1))

#### ROC curves of the features with the highest performance
selectedIndices <- rfCorrRFE$pred$Variables == 17
trainROC <- plot.roc(predictor = rfCorrRFE$pred$Cancer[selectedIndices],
                    x = rfCorrRFE$pred$obs[selectedIndices], legacy.axes = TRUE)

trainPred <- prediction(predictions = rfCorrRFE$pred$Normal[selectedIndices], labels = rfCorrRFE$pred$obs[selectedIndices], label.ordering = rev(levels(c("Cancer", "Normal"))))
trainPerf <- performance(trainPred, "tpr", "fpr")
plot(trainPerf, main = "Title", col = "#1B9E77",lwd = 2, cex.axis = 10)
confusionMatrix(data = rfCorrRFE$pred$pred[selectedIndices],reference = rfCorrRFE$pred$obs[selectedIndices])

#### Model prediction and performance on the testData set
p <- predict(rfCorrRFE, testData)
p.prob <- predict(rfCorrRFE, testData, type = "prob")
confusionMatrix(data = p$pred, reference = testData$class)
testROC <- plot.roc(predictor = p.prob$Cancer,
                   x = testData$class, legacy.axes = TRUE)
testROC
testPred <- prediction(predictions = p.prob$Normal, labels = testData$class, label.ordering = rev(levels(c("Cancer", "Normal"))))
testPerf <- performance(testPred, "tpr", "fpr")
plot(testPerf, main="Title", col ="#7570B3",lwd=2, cex.axis=10)

### ROC curve plot of training and testing results of the model 
plot.ROC(trainModel = rfCorrRFE, testModel = p, testDataClass = testData$class, Title = "ROC Curve - Random Forest")

#### Model testing on an independent cohort of isolated monocytes from Periodontitis
# load the periodontitis dataset
periodontitisTestData = periodontitis.data(predictors = rfCorrRFE$optVariables)
# predict using the trained model
pr <- predict(object = rfCorrRFE, newdata = periodontitisTestData)
confusionMatrix(data = pr$pred, reference = periodontitisTestData$class)



