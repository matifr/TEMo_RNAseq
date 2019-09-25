##------------------------------------------
## Author: Matina Fragkogianni
## Date: 7-9-2018
##
## The main script that performs differential expression analysis (DEA) of 
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
library(ggrepel)
library(biomaRt)
library(pheatmap)
library(sva)

##---------------------------------------------------
# Data import and differential expression analysis
##---------------------------------------------------

pData <- read.delim(file = "phenoData/ALL_pData.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
pData <- pData[which(pData$Batch_info %in% c("6", "4", "5")),]
pData <- pData[which(pData$Condition %in% c("Breast_cancer","Normal", "Endometrial_cancer")),]

#### Read data and build the DGE object
d <- readDGE(pData$Filenames, path = "Counts/Batch 1_3_4_5_6/", columns = c(1,2), pData$Condition)

d$samples$batch <- pData$Batch_info
d$samples$protocol <- pData$Sorting.protocol
d$samples$subtype <- pData$Type

#### Filtering of lowly expresse genes and recalculation of library size
keep <- rowSums(cpm(d) > 1) >= 35 
d <- d[keep,]
colnames(d) <- pData$Sample_name
d$samples$lib.size <- colSums(d$counts)

#### Normalization and log2 tranformation
d <- calcNormFactors(d, method = "upperquartile")
lcpm <- cpm(d,log=TRUE,prior.count = 5)

##-------------------
# PCA plot
##-------------------

df_pca <- prcomp(t(lcpm))
df_out <- as.data.frame(df_pca$x)
df_out$sampleConditions <- d$samples$group
df_out$batch <- as.character(d$samples$batch)
df_out$protocol <- as.character(d$samples$protocol)

theme <- theme(panel.background = element_blank(),panel.border = element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background = element_blank(),axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour = "black"), axis.tex = element_text(size=12), axis.ticks = element_line(colour = "black"),
             plot.margin = unit(c(1,1,1,1),"line"), legend.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 12, face = "bold")
             ,axis.title.y = element_text(color="black", size = 12, face = "bold")) 
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p <- ggplot(df_out,aes(x=PC1,y=PC2,color=sampleConditions, shape = sampleConditions))
p <- p + geom_point(size = 3, stroke = 0.8) + scale_shape_manual(values = c(0,1,2)) + theme + xlab(percentage[1]) + ylab(percentage[2]) + geom_text_repel(aes(label = mcolnames(lcpm))) 
p + scale_color_manual(values = c("#D95F02", "#8B008B", "#1B9E77", "#FF1493" )) # for cancer types pca
p + scale_color_manual(values = c("#22FFB3", "#6EA0FF", "#FFBF44")) #batch

#### Batch effect correction using comBat
combatEdata <- run.ComBat(dge = d, mat = lcpm)

##---------------------------------------------------
# Heatmap of top 100 most variable genes after batch effect correction
##---------------------------------------------------
topVarGenes <- head(order(rowVars(combatEdata), decreasing = TRUE),100)
matrix <- combatEdata[topVarGenes,]

df <- data.frame(Condition = pheno$Group)
rownames(df) <- colnames(combatEdata)
names(df)<- c("Condition")

Condition        <- c("#D95F02", "#8B008B", "#1B9E77")
names(Condition) <- levels(pheno$Group)

anno_colors <- list(Condition = Condition)
pheatmap(matrix, annotation_col = df, annotation_colors = anno_colors,  scale = "row",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         show_rownames= FALSE, show_colnames = F)

##---------------------------------------------------
# Differential expression analysis (DEA)
##---------------------------------------------------

batch <- d$samples$batch
group <- factor(d$samples$group)
group <- relevel(group, ref="Normal")
design <- model.matrix(~0 +group)

cont.matrix <- makeContrasts(
  BCvsNorm <- groupBreast_cancer-groupNormal,
  EndvsNorm <- groupEndometrial_cancer - groupNormal,
  BCvsEndo = groupBreast_cancer-groupEndometrial_cancer,
  levels=design)

#### Fit the linear model for differential expression analysis
vfit <- lmFit(combatEdata, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
efit <- eBayes(vfit, robust = TRUE)
plotSA(efit, main = "Final model: Mean−variance trend")

#### Filter the gene table to identify significant genes
tt <- topTable(efit, adjust.method ='BH', coef = 1, number=999999)
sigGenes <- tt[tt$adj.P.Val <= 0.05,] 
sigGenes_up_down <- tt[tt$adj.P.Val <= 0.05 & (tt$logFC >= 1.5  | tt$logFC <= -1.5),]
up.regulated.genes <- sigGenes[sigGenes$logFC >= 1.5, ]
down.regulated.genes <- sigGenes[sigGenes$logFC <= -1.5, ]

dim(sigGenes_up_down)
dim(up.regulated.genes)
dim(down.regulated.genes) 

##---------------------------------------------------
# Barplot of most interesting genes expressed by monocytes
##---------------------------------------------------
gene_list <- read.delim("Data/gene_list_monocytes.txt", stringsAsFactors = FALSE)
gene_list_2 <- merge(as.data.frame(test_vlookup), gene_list, by.x="external_gene_name", by.y="Gene_name") 

# Sort by Group
gene_list_2 <- with(gene_list_2, gene_list_2[order(Group, external_gene_name, logFC),]) 
gene_list_2$external_gene_name <- factor(gene_list_2$external_gene_name, levels=gene_list_2$external_gene_name)
gene_list_2 <- gene_list_2[gene_list_2$logFC >=1 | gene_list_2$logFC <= -1,]

require(scales)
n <- length(levels(gene_list_2$Group)) # number of colors
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 1, direction = 1)(n)[c(1, 10, 5)] # color palette in random order

my.cols <- brewer.pal(5, "Set1")
my.cols2 <- brewer.pal(5, "Set1")
my.cols[1:4] <- c("#d95f02","#e6ab02","#377eb8","#984ea3")
my.cols2[1:4] <- c("#fb8072","#bebada","#80b1d3","#b3de69")

ggplot(gene_list_2, aes(x = external_gene_name, y = logFC, fill = factor(Group))) + theme_classic() + labs(fill = "") +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") + scale_fill_manual(values = my.cols2) +
  xlab("") +
  ylab("LogFC") + ylim(c(-5,5)) + 
  theme(axis.text=element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(text = element_text(size=16)) +
  theme(axis.line = element_line(size = 0.5, colour = "black"))
