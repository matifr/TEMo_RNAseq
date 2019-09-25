##------------------------------------------
## Author: Matina Fragkogianni
## Date: 1-10-2018
##
## The main script that performs differential expression analysis (DEA) of non-classical
## monocytes from human breast and endometrial cancer samples
## 
##------------------------------------------

##------------------------------------------
# Load libraries
##------------------------------------------
library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)
require(ggrepel)
library(biomaRt)
library(pheatmap)
library(GGally)
library(genefilter)

##------------------------------------------
# Source files
##------------------------------------------
source(file = "Classification_Analysis/tidy_code/run.ComBat.R")

##------------------------------------------
# Data import and normalization 
##------------------------------------------
pData <- read.delim(file = "phenoData/ALL_pData.txt", header = TRUE, stringsAsFactors = FALSE)
pData <- pData[which(pData$Monocyte_population == "Non_Classical"),]
pData <- pData[-which(pData$Sample_names == "mon6"),]

#### Read data and build the DGE object
d <- readDGE(pData$Filenames, path = "Counts/Batch 1_3_4_5_6//", columns = c(1,2), pData$Condition)
d$samples$Condition <- as.factor(pData$Monocyte_population)
d$samples$batch <- pData$Batch_info

#### Filtering of lowly expresse genes and recalculation of library size
keep <- rowSums(cpm(d) > 1) >= 5   
d <- d[keep,]
colnames(d) <- pData$Sample_names
d$samples$lib.size <- colSums(d$counts)

#### Normalization and log2 tranformation
d <- calcNormFactors(d, method = "upperquartile")
lcpm <- cpm(d,log=TRUE,prior.count=5)


#### Batch effect correction using comBat
combatEdata <- run.ComBat(dge = d, mat = lcpm)


##-------------------
# PCA plot
##-------------------

df_pca <- prcomp(t(combatEdata))
df_out <- as.data.frame(df_pca$x)
df_out$pop <- d$samples$Condition
df_out$sampleConditions <- d$samples$group
df_out$batch <- d$samples$batch

theme <- theme(panel.background = element_blank(),panel.border = element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background = element_blank(),axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour = "black"), axis.text = element_text(size = 12), axis.ticks = element_line(colour = "black"),
             plot.margin = unit(c(1,1,1,1),"line"), legend.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 12, face = "bold")
             ,axis.title.y = element_text(color = "black", size = 12, face = "bold")) 
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep = "") )

p <- ggplot(df_out,aes(x = PC1,y = PC2,color = sampleConditions, shape = pop))
p <- p+geom_point(size = 3, stroke = 0.8)  + scale_shape_manual(values = c(18,2)) + theme + xlab(percentage[1]) + ylab(percentage[2]) 
p + scale_color_manual(values = c( "#E41A1C", "#984EA3", "#377EB8"))

##---------------------------------------------------
# Differential expression analysis (DEA)
##---------------------------------------------------

batch <- pData$Batch_info
group <- factor(d$samples$group)
group <- relevel(group, ref = "Normal")
design <- model.matrix(~0 + group)


cont.matrix <- makeContrasts(
  BCvsNorm = groupBreast_cancer-groupNormal,
  ENDvsNorm = groupEndometrial_cancer-groupNormal,
  BCvsEND = groupBreast_cancer-groupEndometrial_cancer,
  levels=design)

vfit <- lmFit(combatEdata, design)
vfit <- contrasts.fit(vfit, contrasts = cont.matrix)
efit <- eBayes(vfit, robust = TRUE)
plotSA(efit, main="Final model: Mean−variance trend")

tt <- topTable(efit, adjust.method ='BH', coef = 1,number = 999999)
sigGenes <- tt[tt$adj.P.Val <= 0.05,] #get seperataly the significant genes
sigGenes_up_down <- tt[tt$adj.P.Val <= 0.05 & (tt$logFC >=1.5  | tt$logFC <= -1.5),] #get seperataly the significant genes
up.regulated.genes <- sigGenes[sigGenes$logFC >= 1.5, ]
down.regulated.genes <- sigGenes[sigGenes$logFC <= -1.5, ]

dim(sigGenes_up_down)
dim(up.regulated.genes)
dim(down.regulated.genes)

##---------------------------------------------------
# Heatmap of significantly differentially expressed genes
##---------------------------------------------------

matrix <- combatEdata[which(rownames(combatEdata) %in% rownames(sigGenes_up_down)),]
df <- data.frame(Condition = d$samples$group)
rownames(df) <- colnames(combatEdata)
names(df) <- c("Condition") 


Condition <- c( "#E41A1C","#984EA3" ,"#377EB8")
names(Condition) <- levels(d$samples$group)

anno_colors <- list(Condition = Condition)
pheatmap(matrix, annotation_col = df, annotation_colors = anno_colors,  scale = "row",
         clustering_distance_cols = "euclidean", clustering_method = "complete",show_rownames = FALSE, show_colnames = FALSE)








