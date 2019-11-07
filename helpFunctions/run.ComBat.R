##------------------------------------------
## Author: Matina Fragkogianni
## Date: 28-9-2018
##
## A function that runs batch effect correction on a gene expression matrix
## @param dge a DGE object
## @param mat a log2cpm normalised gene expression matrix o
## @return: Returns the batch effect corrected matrix of gene expression
## 
##------------------------------------------

run.ComBat <- function(dge,mat){
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(sva)
  
  #### Create the pheno data frame
  pheno <- data.frame(Snames = dge$samples$files, Group = dge$samples$group, Batch = dge$samples$batch)
  
  #### Model matrix including the group but not the sorting protocol
  modcombatCoef <- model.matrix(~as.factor(Group), data=pheno)
  combatEdata <- sva::ComBat(dat = mat, batch = pheno$Batch, mod = modcombatCoef, par.prior = TRUE, prior.plots = FALSE)
  
  return(combatEdata)
}