library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)
require("ggrepel")
library(biomaRt)
library(pheatmap)
library(GGally)
library(genefilter)

pData = read.delim(file = "phenoData/ALL_pData.txt", header = TRUE, stringsAsFactors = FALSE)
pData = pData[which(pData$Monocyte_population == "Non_Classical"),]
pData = pData[-which(pData$Sample_names == "mon6"),]
#write.table(pData, file = "~/Desktop/pData_noncl_mono.txt", sep = "\t", row.names = FALSE)

d =readDGE(pData$Filenames, path = "Counts/Batch 1_3_4_5_6//", columns=c(1,2), pData$Condition)
d$samples$Condition = as.factor(pData$Monocyte_population)
d$samples$batch = pData$Batch_info
table(rowSums(d$counts==0)==19)
table(d$samples$group)

keep = rowSums(cpm(d)>1) >= 5   
d = d[keep,]
table(keep)

colnames(d) = pData$Sample_names
d$samples$lib.size = colSums(d$counts) #recalculate the library size after filtering
d <- calcNormFactors(d, method = "upperquartile")

lcpm = cpm(d,log=TRUE,prior.count=5)

pheno = data.frame(Snames = colnames(d), Group = d$samples$group, Batch = d$samples$batch)

modcombat_coef = model.matrix(~as.factor(Group), data=pheno)
combat_edata = ComBat(dat = lcpm, batch = pheno$Batch, mod = modcombat_coef, par.prior = TRUE, prior.plots = FALSE)

#### Exploratory analysis on monocyte populations of batch 1
df_pca <- prcomp(t(combat_edata))
df_out = as.data.frame(df_pca$x)
df_out$pop = d$samples$Condition
df_out$sampleConditions = d$samples$group
df_out$batch = d$samples$batch

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"), axis.text=element_text(size=12), axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"), legend.text=element_text(size=12), axis.title.x = element_text(color="black", size=12, face="bold")
             ,axis.title.y = element_text(color="black", size=12, face="bold")) 
theme_update(plot.title = element_text(hjust = 0.5, face="bold"))
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=sampleConditions, shape = pop))
p<-p+geom_point(size = 3, stroke = 0.8)  + scale_shape_manual(values = c(18,2)) + theme + xlab(percentage[1]) + ylab(percentage[2]) 
p + scale_color_manual(values=c( "#E41A1C", "#984EA3", "#377EB8"))


# topVarGenes <- head(order(rowVars(combat_edata), decreasing=TRUE),100)
# matrix <- combat_edata[topVarGenes,]
# 
# df <- data.frame(Condition = d$samples$group)
# rownames(df) = colnames(combat_edata)
# names(df) = c("Condition") 
# 
# 
# Condition        <- c( "#E41A1C","#984EA3" ,"#377EB8")
# names(Condition) <- levels(d$samples$group)
# 
# anno_colors <- list(Condition = Condition)
# pheatmap(matrix, annotation_col = df, annotation_colors = anno_colors,  scale = "row",clustering_distance_cols = "euclidean", clustering_method = "complete")


batch = pData$Batch_info
group=factor(d$samples$group)
group = relevel(group, ref="Normal")
design = model.matrix(~0 + group)


cont.matrix = makeContrasts(
  BCvsNorm = groupBreast_cancer-groupNormal,
  ENDvsNorm = groupEndometrial_cancer-groupNormal,
  BCvsEND = groupBreast_cancer-groupEndometrial_cancer,
  levels=design)

#v = voom(d, design = design,plot = TRUE)
#v
vfit <- lmFit(combat_edata, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
efit <- eBayes(vfit, robust = TRUE)
plotSA(efit, main="Final model: Mean−variance trend")

tt = topTable(efit, adjust.method ='BH', coef = 1,number=999999)
sigGenes= tt[tt$adj.P.Val <= 0.05,] #get seperataly the significant genes
sigGenes_up_down= tt[tt$adj.P.Val <= 0.05 & (tt$logFC >=1.5  | tt$logFC <= -1.5),] #get seperataly the significant genes
up.regulated.genes = sigGenes[sigGenes$logFC >= 1.5, ]
down.regulated.genes = sigGenes[sigGenes$logFC <= -1.5, ]

dim(sigGenes_up_down)
dim(up.regulated.genes)
dim(down.regulated.genes)


source(file = "Functions/remove.dots.R")
rownames(sigGenes_up_down) = remove.dots(sigGenes_up_down)
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "grch37.ensembl.org")
all_genes= getBM(filters = "ensembl_gene_id", attributes= c("external_gene_name","ensembl_gene_id"),
                 values = rownames(sigGenes_up_down), mart = ensembl)

test_vlookup=merge(as.data.frame(sigGenes_up_down), all_genes, by.x="row.names", by.y="ensembl_gene_id") 

write.table(test_vlookup, file = "~/Desktop/ENDO_Norm_sig_logFC_genes.txt", sep = "\t",row.names = F)


matrix <- combat_edata[which(rownames(combat_edata) %in% rownames(sigGenes_up_down)),]
df <- data.frame(Condition = d$samples$group)
rownames(df) = colnames(combat_edata)
names(df) = c("Condition") 


Condition        <- c( "#E41A1C","#984EA3" ,"#377EB8")
names(Condition) <- levels(d$samples$group)

anno_colors <- list(Condition = Condition)
pheatmap(matrix, annotation_col = df, annotation_colors = anno_colors,  scale = "row",
         clustering_distance_cols = "euclidean", clustering_method = "complete",show_rownames= FALSE, show_colnames = FALSE)








