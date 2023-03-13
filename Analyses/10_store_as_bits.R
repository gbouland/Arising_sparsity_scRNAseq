library(reshape2)
library(patchwork)
library(ggplot2)

getReduction <- function(counts){
  library(bit)
  library(Matrix)
  library(magrittr)
  library(Seurat)
  counts <- counts
  counts <- counts[,colSums(counts)>0]
  colnames(counts) <- sprintf("colNo_%s",1:ncol(counts))
  rownames(counts) <- sprintf("rowNo_%s",1:nrow(counts))
  seur <- CreateSeuratObject(counts = counts)
  seur <- SCTransform(seur)
  if(ncol(counts)>25000){    
    splits <- split(1:ncol(counts), ceiling(1:ncol(counts)/25000))
    newMatr <- matrix(nrow = nrow(counts),ncol = ncol(counts))
    for(x in 1:length(splits)){
      message(sprintf("%s of %s",x,length(splits)))
      newMatr[,splits[[x]]] <- as.matrix(counts[,splits[[x]]])
    }
    counts <- newMatr 
    rm(newMatr)
  }else{
    counts <- as.matrix(counts)
  }
  binary <- (counts >=1)
  List <- vector("list", ncol(binary))
  for(i in 1:length(List)){
    message(i)
    newBit <- bit(nrow(binary))
    newBit[which(binary[,i])] <- T
    List[[i]] <- newBit
  }
  bitStored <- format(object.size(List), units = "Mb")%>% gsub(pattern = " Mb",replacement = "") %>% as.numeric()
  scTransform <- format(object.size(seur@assays$SCT), units = "Mb")%>% gsub(pattern = " Mb",replacement = "") %>% as.numeric()
  toSave <- data.frame("bits" = bitStored, "countsNormalized" = scTransform)    
  return(toSave)
}
getReduction(counts)## For every dataset this was added to the full summary

summary_data <- read.csv2("./Data/fullSummary.csv")
plotSizes <- summary_data[,c("dataset","bitSize","lognormSize","scTransSize")]
plotData <- melt(plotSizes,id.vars = "dataset")
plotData$dataset <- factor(plotData$dataset, levels = summary_data[order(summary_data$scTransSize),"dataset"])
a <- ggplot(plotData,aes(value,dataset, fill = variable)) + geom_bar(stat = "identity",position = position_dodge()) + theme_minimal()
plotSizes$bit_sc <- plotSizes$scTransSize / plotSizes$bitSize
plotSizes$bit_ln <- plotSizes$lognormSize / plotSizes$bitSize
plotSizes2 <- plotSizes[,c("dataset","bit_sc","bit_ln")]
plotData2 <- melt(plotSizes2,id.vars = "dataset")
plotData2$dataset <- factor(plotData2$dataset, levels = plotSizes2[order(plotSizes2$bit_sc),"dataset"])
b <- ggplot(plotData2,aes(value,dataset, fill = variable)) + geom_bar(stat = "identity",position = position_dodge()) + theme_minimal() + viridis::scale_fill_viridis(discrete = T)

a | b


