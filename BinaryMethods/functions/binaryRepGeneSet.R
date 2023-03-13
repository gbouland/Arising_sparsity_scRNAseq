##Get a representative set of genes
##binary    binarized single-cell data
##ngenes    number of genes
binaryRepGeneSet <- function(binary, ngenes = 2000){
  ranges <- seq(0.1,0.9,0.05)
  proportions <- round(-log(ranges[1:length(ranges)-1]) / sum(-log(ranges[1:length(ranges)-1])) * ngenes)
  dr <- rowMeans(binary)
  sel_genes <- lapply(1:(length(ranges)-1), function(i){
    genes <- names(dr[dr>= ranges[i] & dr<= ranges[i+1]])
    sample(genes,min(c(proportions[i],length(genes))))
  }) |> unlist()
  return(sel_genes)
}