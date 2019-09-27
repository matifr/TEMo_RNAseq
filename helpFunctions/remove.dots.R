
remove.dots = function(data){
  g = rownames(data)
  g = strsplit(g,"[.]")
  genes_names = as.character(length(g))
  
  for (i in 1:length(g)){
    genes_names[i] = g[[i]][1]
  }
  
  return (genes_names)
}