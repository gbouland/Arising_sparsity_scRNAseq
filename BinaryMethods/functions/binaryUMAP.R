## UMAP based on binarized single-cell data
## binary     binarized single-cell data
## reduction  dimensionality reduction method (PCA, BFA or JEV)
## nc         number of components to use
binaryUMAP <- function(binary, reduction, nc = 10, n_neighbors = 10, min_dist = 0.3, spread = 1, metric = "euclidean"){
  require(uwot)
  comps <- switch(reduction,
         "PCA" = {
           binaryPCA(binary,nc)
         },
         "BFA" = {
           binaryBFA(binary,nc)
         },
         "JEV" = {
           JaccardEV(binary,nc)
         })
  UMAP <- uwot::umap(comps ,n_neighbors = n_neighbors, min_dist= min_dist,spread = spread,metric = metric)
  return(UMAP)
}