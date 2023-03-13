##Calculating eigenvectors from jaccard similarity matrix
## binary   is binarized single cell matrix
## n        is number of components
JaccardEV <- function(binary, n=10){
  require(locStra)
  require(Matrix)
  require(svd)
  sparseM <- Matrix(binary,sparse=TRUE)
  jac <- locStra::jaccardMatrix(sparseM, useCpp = TRUE, sparse = TRUE)
  jac_eigen <- svd::trlan.eigen(jac,neig = n)
  jac_eigen <- jac_eigen$u
  rownames(jac_eigen) <- colnames(binary)
  colnames(jac_eigen) <- paste0("Jaccard_",1:n)
  return(jac_eigen) 
}