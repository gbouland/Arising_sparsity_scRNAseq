## Recover magnitude of expression from binarized single-cell RNAseq dataset
## binary     binarized single cell data
## n          number of neighbours to use for recovery
recoverMagnitude <- function(binary, n = 2){
  require(locStra)
  require(Matrix)
  require(matrixStats)
  jac <- jaccardMatrix(binary, useCpp = TRUE, sparse = TRUE)
  colnames(jac) <- colnames(binary)
  rownames(jac) <- colnames(binary)
  recovered <- lapply(1:nrow(jac),function(i){
    message(sprintf("%s of %s",i,nrow(jac)))
    row <- jac[i,]
    neighbours <- sort(row,decreasing = T)[1:n]
    sel <- names(neighbours)
    weights <- unname(neighbours) / sum(neighbours)
    new_vec <- rowWeightedMeans(as.matrix(binary[,sel]),w = weights)
    return(new_vec)
  })
  recovered <- recovered |> do.call(what = "cbind")
  colnames(recovered) <- colnames(binary)
  recovered <- Matrix(recovered, sparse = TRUE)
  norm <- t(t(recovered)/colSums(recovered))* 10000
  recovered <- log(norm + 1)
  recovered[binary ==0] <- 0
  rownames(recovered) <- rownames(binary)
  return(recovered)
}
