## Leiden clustering based on binarized single-cell data
## binary     binarized single-cell data
## reduction  dimensionality reduction method (PCA, BFA or JEV)
## nc         number of components to use
## resolution resolution of leiden clustering
binaryClustering <- function(binary, reduction, nc = 10, resolution = 0.5){
  require(igraph)
  require(dbscan)
  require(locStra)
  require(leiden)
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
  
  
  distance <- dist(comps,method = "euclidean")
  SNN <- dbscan::sNN(distance, k=30)
  IDs <- colnames(binary)
  sNN_graph <- lapply(rownames(SNN$id),function(i){
    neighbours <- IDs[SNN$id[i,]]
    out <- data.frame(source = i,target = neighbours)
  }) %>% do.call(what = rbind)
  links <- sNN_graph
  nodes <- data.frame(nodes = unique(links$source))
  graph <- igraph::graph_from_data_frame(d=links, vertices=nodes, directed=F)
  partition <- leidenAlg::find_partition(graph,edge_weights = rep(1,nrow(links)),resolution = resolution)
  out <- data.frame("ID" = colnames(binary),
                    "cluster" = partition)
  return(out)
}