##------------------------------------------
## Author: Matina Fragkogianni
## Purpose of script: Translate Ensembl IDs to Gene symbols using the biomart database 
## Input: String vector (subset)
## Output: Prints ensembl IDs and corresponding gene symbols on the stdout
## Date: 28-9-2018
##------------------------------------------

ensembl.to.GS <- function(subset){

  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(biomaRt)
  
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "grch37.ensembl.org")
  all_genes= getBM(filters = "ensembl_gene_id", attributes= c("external_gene_name","ensembl_gene_id", "description"),
                   values = subset, mart = ensembl)
  
  print(all_genes)
  
}