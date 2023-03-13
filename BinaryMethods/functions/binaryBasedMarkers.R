## Find markers based on binarized single-cell data
##binary    binarized single-cell data
##clusters  clusters for which markers have to be found
findMarkers <- function(binary, clusters){
  cluster_matrix <- lapply(unique(clusters),function(x){
    as.integer(x == clusters)
  }) |> do.call(what = "cbind")
  colnames(cluster_matrix) <- unique(clusters)
  rownames(cluster_matrix) <- colnames(binary)
  res <- cor(x = t(as.matrix(binary)), y = as.matrix(cluster_matrix))
  t_matr <- res / sqrt((1-res^2)/(nrow(cluster_matrix) -2))
  p_matr <- 2*pt(t_matr, nrow(cluster_matrix)-2, lower=FALSE)
  p_matr[p_matr >1] <-1
  markers <- lapply(colnames(res),function(x){
    tmpRes <- res[order(res[,x],decreasing = T),]
    tmpP <- p_matr[order(res[,x],decreasing = T),]
    data.frame("gene" = rownames(tmpRes),
               "cluster" = x,
               "cor" = tmpRes[,x],
               "p"  = tmpP[,x])
  })
  names(markers) <- colnames(res)
  return(markers)
}