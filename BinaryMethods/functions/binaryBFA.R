##Wrapper around scBFA (binary factor analysis)
## binary   is binarized single cell matrix
## n        is number of components
binaryBFA <- function(binary, n=10){
  sce <- SingleCellExperiment(assay = list(logcounts = binary))
  bfa_model = scBFA(scData = sce, numFactors = n)
  cell_loadings <- bfa_model$ZZ
  rownames(cell_loadings) <- colnames(binary)
  colnames(cell_loadings) <- paste0("scBFA_",1:n)
  return(cell_loadings)
}







