##------------------------------------------
## Author: Matina Fragkogianni
## Date: 7-11-2018
##
## This is the main script to perform permutation testing and plot a histogram 
## 
##------------------------------------------

##------------------------------------------
# Source files
##------------------------------------------
source(file = "helpFunctions//random.class.perm.test.R")

#### Run permutation testing for the training data 
random.classifiers <- random.class.perm.test(data_fr, rfCorrRFE$optVariables, 1000)

#### get summary statistics for the accuracy of the random classifiers
summary(random.classifiers$TrainAccuracy)
sd(random.classifiers$TrainAccuracy)
quantile(random.classifiers$TrainAccuracy, probs = c(.025, 0.975))

#### get significance of the observed accuracy compared to random classifiers 
(sum( c(random.classifiers$TrainAccuracy, 0.8509) >= 0.8509) +1)/ (1000 +1) # one tailed test for training

# print density plot and histogram for training accuracy
df <- data.frame("Performance" = random.classifiers$TrainAccuracy)

ggplot(df, aes(x = Performance, color = Performance)) +
  geom_histogram(color = "#1B9E77", fill = "#1B9E77",bins = 20,alpha = 0.5,position = "identity") +
  geom_vline(data = df, aes(xintercept=0.8509, color = "red"),linetype = "solid") +
  labs(title = "Random signatures (training set)", y = "Histogram") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face = "bold", size = 20)) 
 
