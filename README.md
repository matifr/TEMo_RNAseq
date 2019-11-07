# Transcriptomic analysis of human monocytes for breast cancer diagnosis

## Table of contents
* [Background and motivation ](#Background-and-motivation )
* [Project Aims ](#Project-Aims )
* [How to run the analysis](#How-to-run-the-analysis)
	* [Differential expression analysis](#Differential-expression-analysis)
	* [Classification](#Classification)

## Background and motivation 

Breast and endometrial cancer are two of the most commonly diagnosed cancers in women in the UK. There is increasing evidence that early detection of tumours improves survival rates. Mammography is the most accurate and non-invasive method to date for early detection of breast cancer, however its efficacy drops in women with higher breast density, resulting in overdiagnosis and potentially unnecessary interventions (Kolb et al. 2002). Currently there is no available screening method for early detection of endometrial cancer. Therefore, there is an urgent need for a non-invasive method for early detection of breast and endometrial cancer. 
The advent of high-throughput microarray technologies has resulted in an increased interest in identification of multigene signatures for diagnosis and prognosis of cancer.Monocytes are circulating blood cells that migrate to tissues and give rise to macrophages. Recent studies have revealed a pro-tumoral profile for these cells in renal carcinoma and colorectal cancer (Chittezhath et al. 2014; Hamm et al. 2015). Consequently, this study investigated the transcriptomes of circulating monocytes in breast and endometrial cancer in order to develop a blood-based, non-invasive diagnostic test for the detection of cancer. This test would act as a clinical tool for efficient screening of patients at risk. 

## Methods

### Classification 
Recursive feature elimination (RFE) is a wrapper feature selection method that starts by fitting the classification model to all features and then each feature is ranked for its importance to the model. We assume that S is a sequence of values that indicates the number of features to be retained. At each iteration, the Si top ranked features are used to refit the model and the subset of genes with the highest accuracy is selected. I used the RFE feature selection with resampling and a Random forest (RF) classifier as implemented in the caret package in Bioconductor (Kuhn, 2015). In brief, RF is an ensemble learning method that constructs many decision trees (forest) and uses a majority vote to make predictions (Breiman, 2001). Assuming N number of samples in the training set, the algorithm creates a subset of the data of the N samples with replacement, and at each node, m number of genes are selected at random from all the genes M. The variable m that gives the best split is selected and used for a binary split. This procedure is repeated for each node until the tree is grown to terminal nodes k.

Before model training and feature selection, the breast TEMo dataset of n = 77 samples (32 healthy, 45 breast cancer) was filtered for lowly expressed genes (CPM > 1 across conditions), normalized using upper-quantile normalization, log2 transformed using the cpm() function and corrected for batch effects using ComBat. Then, the dataset was split into training (70%, n = 55, 32 healthy, 23 cancer samples) and testing (30%, n = 22, 13 healthy, 9 cancer samples). While the test set was kept on the side, the training set was used for RFE-RF model training and feature selection using 5 times repeated 10-fold cross-validation (CV). In short, the training samples were randomly partitioned into k (k = 10) subsamples. Out of the k sub-samples, one is kept for testing the classifier, and the remaining k-1 subsamples are used for training the classifier. Then the whole process is repeated n= 5 times. For the training set, overall accuracy, sensitivity and specificity were calculated averaging the CV predictions for the optimal subset. The model with the highest average accuracy was selected as optimal. The optimal model was the fitted on the test set that had not been used for training or feature selection. For the RFE-RF model, subsets of features were selcted ranging between 2 to 30, variable ntree was set to 500 trees and variable mtry was calculated as √p, where p is the number of features used during training the model. The receiver-operating characteristic (ROC) curves were drawn using the ROCR package in R (Sing et al., 2005). To determine the accuracy rates of the classifiers that can be obtained by chance, a RF model using the 17 selected genes was trained with permuted class labels. This process was performed during training 1,000 times using 5 x 10-fold CV. The accuracy of the 1,000 random classifiers was recorded. The p value was calculated by counting the accuracy of the random classifiers that achieved similar or higher total accuracy compared to the observed accuracy of the RFE-RF model on the training data.

## Project Aims 
1. To compare the transcriptional profiles of circulating monocytes coming from breast and endometrial cancer patients to healthy individuals as well as between cancer types.

2. To identify a diagnostic signature and develop a classification model for stratification of cancer patients.

3. To validate the robustness of the diagnostic signature on independent patient cohorts as well as datasets of monocytes from different cancer types.

## How to run the analysis

### Differential expression analysis
To perform differential expression analysis run <code> TEMo_DEA.R </code>. This script loads the raw counts, performs *upperquantile* normalization, batch effect correction and differential expression using the *limma* package. It produces a PCA plot, a heatmap of the top 100 most variable genes between conditions as well as a barplot with a selection of interesting genes.

### Classification
coming soon..
