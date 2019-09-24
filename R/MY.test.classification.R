# Random classifiers
source(file = "Classification_Analysis/tidy_code/random.class.perm.test.R")
random.classifiers = random.class.perm.test(data_fr, rf.corr.rfe$optVariables, 1000)
#save(random.classifiers, file="Classification_Analysis/tidy_code/random.classifiers_1000_permutations.RData")
load("~/repos_not_updated/nature_paper_code/Monocytes_analysis_2018/Classification_Analysis/tidy_code/random.classifiers_1000_permutations.RData")
summary(random.classifiers$TrainAccuracy)
sd(random.classifiers$TrainAccuracy)
quantile(random.classifiers$TrainAccuracy, probs = c(.025, 0.975))
quantile(random.classifiers$TestAccuracy, probs = c(.025, 0.975))
summary(random.classifiers$TestAccuracy)

#histogram and p value Testing
#hist(random.classifiers$TestAccuracy)
#(sum(random.classifiers$TestAccuracy > 0.8182) + 1) / (1000 + 1)  # one tailed test for testing
#(sum( c(random.classifiers$TestAccuracy, 0.8182) >= 0.8182) +1)/ (1000 +1) # one tailed test for training

# p value training
(sum( c(random.classifiers$TrainAccuracy, 0.8509) >= 0.8509) +1)/ (1000 +1) # one tailed test for training

# print density plot and histogram for trainnig accuracy
df = data.frame("Performance" = random.classifiers$TrainAccuracy)
p <- ggplot(df, aes(x=Performance)) + 
  geom_density(fill="#56B4E9", alpha=0.4) 
p+ theme(legend.position="none") + 
  geom_vline(data=df, aes(xintercept=0.8182, color="red"),linetype="solid") + xlim(-0.5,1.3) +
  labs(title="Random signatures", y = "Density") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold", size=20)) 

ggplot(df, aes(x=Performance, color = Performance)) +
  geom_histogram(color="#1B9E77", fill="#1B9E77",bins = 20,alpha=0.5,position="identity") +
  geom_vline(data=df, aes(xintercept=0.8509, color="red"),linetype="solid") +
  labs(title="Random signatures (training set)", y = "Histogram") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold", size=20)) 


# print density plot and histogram for test accuracy
df = data.frame("Performance" = random.classifiers$TestAccuracy)
p <- ggplot(df, aes(x=Performance)) + 
  geom_density(fill="#56B4E9", alpha=0.4) 
p+ theme(legend.position="none") + 
  geom_vline(data=df, aes(xintercept=0.8182, color="red"),linetype="solid") + xlim(-0.5,1.3) +
  labs(title="Random signatures", y = "Density") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold", size=20)) 
ggplot(df, aes(x=Performance, color = Performance)) +
  geom_histogram(color="#1B9E77", fill="#1B9E77",bins = 20,alpha=0.5,position="identity") + 
  geom_vline(data=df, aes(xintercept=0.8182, color="red"),linetype="solid") +
  labs(title="Random signatures (testing set)", y = "Histogram") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold", size=20)) 



########################################################################
########################################################################

# Reproducibilty 

RepTrainAccuracy = c()
RepTestAccuracy =c()
data_rep = data_fr[,rf.corr.rfe$optVariables]
data_rep$class = data_fr$class
source(file = "Classification_Analysis/tidy_code/run.RF.R")
for(i in 1:100){
  
  set.seed(i)
  # Step 1: Get row numbers for the training data
  trainRowNumbers <- createDataPartition(data_rep$class, p = 0.7, list=FALSE)
  
  # Step 2: Create the training  dataset
  trainData <- data_rep[trainRowNumbers,]
  dim(trainData)
  
  
  # Step 3: Create the test dataset
  testData <- data_rep[-trainRowNumbers,]
  dim(testData)
  
  cat("Running RF")
  cat(i)
  cat("\n")
  rf.model = run.RF(trainData = trainData)
  
  train.acc = getTrainPerf(rf.model)[[1]]
  
  RepTrainAccuracy[i] = train.acc
  
  p = predict(rf.model, testData)
  pred.conf = confusionMatrix(data = p, reference = testData$class)
  
  RepTestAccuracy[i] = pred.conf$overall[[1]]
}

summary(RepTrainAccuracy)
summary(RepTestAccuracy)

########################################################################
########################################################################
########################################################################

# 1,000 random signatures of the same size
set.seed(123)
x <- replicate(17, {
  index <- sample(1:ncol(trainData), 1000, replace = FALSE)
})

perf  <- vector("list", 1000)
source(file = "Classification_Analysis/tidy_code/run.RF.R")
for(i in 1:nrow(x)){
  set.seed(123)
  random_data = trainData[,x[i,]]
  random_data$class = trainData$class
  RF.model = run.RF(trainData = random_data)
  
  perf[[i]] = getTrainPerf(RF.model)
}


df = data.frame("Performance" = as.matrix(random.classifiers$TrainAccuracy))#, Mode = c("0.58"), sigPerf = c("0.85"))

p <- ggplot(df, aes(x=Performance)) + 
  geom_density(fill="#56B4E9", alpha=0.4, bw = 0.007) 
p + theme(legend.position="none") + 
  geom_vline(data=df, aes(xintercept=0.8588, color="red"),linetype="solid") + xlim(0.5,1) +
  #geom_vline(data=df, aes(xintercept=Mode, color="red"),linetype="dashed") + 
  labs(title="Random signatures",x="Performance", y = "Density") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold", size=20)) 

TrainAccuracy = c()
for(i in 1:1000){
  TrainAccuracy[i] = perf17_genes[[i]]$TrainAccuracy
}
mean(TrainAccuracy)
max(TrainAccuracy)
min(TrainAccuracy)
densityplot(TrainAccuracy)
quantile(TrainAccuracy, probs = c(.05, 0.1, 0.25, 0.5, 0.75, 0.9 , 0.95))
perf17_genes = perf
#save(perf17_genes, file = "~/Desktop/perf17_genes.RData")
(sum( c(TrainAccuracy, 0.9) >= 0.9) +1)/ (1000 +1) # one tailed test for training

