##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Run batch effect correction
## Input: DGE object (dge), Normalized matrix of gene expression (mat)
## Output: Returns the batch effect corrected matrix of gene expression
## Date: 28-9-2018
##------------------------------------------

run.ComBat <- function(dge,mat){
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(sva)
  
  #### Create the pheno data frame
  pheno = data.frame(Snames = colnames(dge), Group = dge$samples$group, Batch = dge$samples$batch)
  
  #### Model matrix including the group but not the sorting protocol
  modcombatCoef = model.matrix(~as.factor(Group), data=pheno)
  combatEdata = ComBat(dat = mat, batch = pheno$Batch, mod = modcombatCoef, par.prior = TRUE, prior.plots = FALSE)
  
  return(combatEdata)
}