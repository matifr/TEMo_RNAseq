##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Principle component analysis plot for gene expressin data
## Input: Normalized matrix of gene expression (mat), DGE object (dge), smaples colored by condition or batch (option)
## Output: PCA plot
## Date: 28-9-2018
##------------------------------------------
my.pca <- function(mat, dge, option = "condition"){
  
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(stats)
  require(ggplot2)
  require(ggrepel)
  
  df_pca <- prcomp(t(mat))
  df_out = as.data.frame(df_pca$x)
  df_out$sampleConditions = dge$samples$group
  df_out$batch = as.character(dge$samples$batch)
  
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"), axis.text=element_text(size=12), axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"), legend.text=element_text(size=12), axis.title.x = element_text(color="black", size=12, face="bold")
               ,axis.title.y = element_text(color="black", size=12, face="bold")) 
  theme_update(plot.title = element_text(hjust = 0.5, face="bold"))
  percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
  percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  
  if(option %in% "condition"){
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=sampleConditions, shape = batch))
    p<-p+geom_point(size = 3, stroke = 0.8) + scale_shape_manual(values = c(0,1,2))+ theme + xlab(percentage[1]) + ylab(percentage[2]) 
    p + scale_color_manual(values=c("#D95F02", "#1B9E77" ))
  }else{
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=sampleConditions, shape = sampleConditions))
    p<-p+geom_point(size = 3, stroke = 0.8) + scale_shape_manual(values = c(0,1,2))+ theme + xlab(percentage[1]) + ylab(percentage[2]) 
    p + scale_color_manual(values=c("#22FFB3", "#6EA0FF", "#FFBF44")) 
  }
}