## Binary based aggregation on detection rate
## binary   binary data
## by       value by which the binary data need to be aggregated
binaryBasedAggregation <- function(binary, by){
  aggregatedMatrix <- matrix(nrow = nrow(binary),
                             ncol = length(unique(by)))
  newColumns <- unique(by)
  for(i in 1:length(newColumns)){
    aggregatedMatrix[,i] <- rowMeans(binary[,which(by == newColumns[i])])
    
  }
  
  colnames(aggregatedMatrix) <- newColumns
  rownames(aggregatedMatrix) <- rownames(binary)
  return(aggregatedMatrix)
}