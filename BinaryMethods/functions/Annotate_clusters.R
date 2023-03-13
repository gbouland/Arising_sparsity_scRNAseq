annotate_clusters <- function(clusterMarkers, celltypemarkers, ntop = 20){
  clusterMarkers <- lapply(clusterMarkers,function(x){
    return(x$gene[1:ntop])
  }) |> do.call(what = "rbind") |> melt()
  clusterMarkers <- clusterMarkers[,c(1,3)]
  one_hot_cluster <- lapply(unique(clusterMarkers$Var1),function(clusID){
    as.numeric(clusterMarkers$value %in% 
                 clusterMarkers[clusterMarkers$Var1 == clusID,"value"])
  }) |> do.call(what = "cbind")
  rownames(one_hot_cluster) <- clusterMarkers$value
  colnames(one_hot_cluster) <- unique(clusterMarkers$Var1)
  markers <- celltypemarkers[celltypemarkers$markers %in% rownames(one_hot_cluster),]
  oneHotMarkers <- dcast(markers,"markers~cell",fill = 0)
  rownames(oneHotMarkers) <- oneHotMarkers$markers
  oneHotMarkers$markers <- NULL
  oneHotMarkers[oneHotMarkers!=0] <- 1
  oneHotMarkers <- as.matrix(oneHotMarkers)
  one_hot_cluster <- one_hot_cluster[rownames(oneHotMarkers),]
  oneHotMarkers <- apply(oneHotMarkers,2,as.numeric)
  rownames(oneHotMarkers) <- rownames(one_hot_cluster)
  props <- cor(oneHotMarkers,as.matrix(one_hot_cluster))
  props <- t(props)
  celltypes <- apply(props,1,function(x){
    sorted <- sort(x,decreasing = T)
    names(sorted)[1]
  })
  return(celltypes)
}



