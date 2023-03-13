library(reshape2)
library(caret)
library(ggplot2)
AD_set <- readRDS("./Data/prep_AD_set.rds")
AD_samplesheet <- AD_set$samplesheet
table(AD_samplesheet$celltype)
binary <- AD_set$binary
markers <- BRETIGEA::markers_df_human_brain
markers <- markers[markers$markers %in% rownames(binary),]
oneHotMarkers <- dcast(markers,"markers~cell",fill = 0)
rownames(oneHotMarkers) <- oneHotMarkers$markers
oneHotMarkers$markers <- NULL
oneHotMarkers[oneHotMarkers!=0] <- 1
oneHotMarkers <- as.matrix(oneHotMarkers)
binary <- binary[rownames(oneHotMarkers),]
oneHotMarkers <- apply(oneHotMarkers,2,as.numeric)
rownames(oneHotMarkers) <- rownames(binary)
props <- cor(oneHotMarkers,as.matrix(binary))
props <- t(props)
celltypes <- apply(props,1,function(x){
  sorted <- sort(x,decreasing = T)
  names(sorted)[1]
})
AD_samplesheet$binary_cell_types <- celltypes
AD_samplesheet[AD_samplesheet$binary_cell_types == "ast","binary_cell_types"] <- "astro"
AD_samplesheet[AD_samplesheet$binary_cell_types == "end","binary_cell_types"] <- "endo"
AD_samplesheet[AD_samplesheet$binary_cell_types == "mic","binary_cell_types"] <- "mg"
AD_samplesheet[AD_samplesheet$binary_cell_types == "oli","binary_cell_types"] <- "oligo"
AD_samplesheet[AD_samplesheet$binary_cell_types == "opc","binary_cell_types"] <- "OPC"
AD_samplesheet[AD_samplesheet$binary_cell_types == "neu","binary_cell_types"] <- "neuron"
conf <- confusionMatrix(table(AD_samplesheet$celltype,AD_samplesheet$binary_cell_types))
##Median F1-score
median(conf$byClass[,7])
long <- melt(conf$table)
colnames(long) <- c("count-based annotations","binary-based annotations","value")
ggplot(long,aes(`count-based annotations`,`binary-based annotations`, fill = value)) + geom_tile() +
  viridis::scale_fill_viridis() + geom_text(aes(`count-based annotations`,`binary-based annotations`, label = value), color = "white", size = 4) + theme_minimal()




source("./BinaryMethods/functions/binaryUMAP.R")
source("./BinaryMethods/functions/binaryPCA.R")
UMAPs <- binaryUMAP(binary = AD_set$binary, reduction = "PCA", nc = 20)
AD_samplesheet$UMAP1 <- UMAPs[,1]
AD_samplesheet$UMAP2 <- UMAPs[,2]
ggplot(AD_samplesheet,aes(UMAP1,UMAP2,col = celltype)) + geom_point() + theme_minimal()

markersToPlot <- c("oligo" = "PLP1",
                   "astro" = "AQP4",
                   "OPC" = "TNR",
                   "neuron" = "SYNPR",
                   "endo" = "IFI27",
                   "mg" = "PTPRC")
markerExpression <- lapply(markersToPlot,function(x){
  out <- cbind(AD_set$binary[x,],AD_set$normalized[x,])
  colnames(out) <- c(paste0(x,"_binary"),paste0(x,"_normalized"))
  return(out)
}) |> do.call(what = "cbind")

toPlot <- colnames(markerExpression)
AD_samplesheet <- cbind(AD_samplesheet,markerExpression)

plots <- lapply(names(markersToPlot),function(x){
  bin <- toPlot[grepl(markersToPlot[x],toPlot)][1]
  norm <- toPlot[grepl(markersToPlot[x],toPlot)][2]
  a <- ggplot(AD_samplesheet,aes(UMAP1,UMAP2,col = as.factor(get(bin)))) +
    geom_point(size = 0.5,alpha = 0.75) +
    viridis::scale_color_viridis(discrete = T) +
    theme_minimal() + labs(col = "Binarized expression",
                           title = sprintf("%s marker for %s",markersToPlot[x],x))
     
  
  b <- ggplot(AD_samplesheet,aes(UMAP1,UMAP2,col = get(norm))) +
    geom_point(size = 0.5) + viridis::scale_color_viridis() +
    theme_minimal() + labs(col = "Normalized expression" , 
                           title = sprintf("%s marker for %s",markersToPlot[x],x))
  return(a | b)
})



plots[[1]] / plots[[2]] / plots[[3]] / plots[[4]] / plots[[5]] / plots[[6]]






a <- ggplot(AD_samplesheet,aes(BPC1,BPC2, col = celltype)) + geom_point(size = 0.9) + theme_minimal() + labs(x = "BinaryPC1",y = "BinaryPC1") +
  theme(legend.position="top")
b <- ggplot(AD_samplesheet,aes(PC1,PC2, col = celltype)) + geom_point(size = 0.9) + theme_minimal() + labs(x = "PC1",y = "PC2") +
   theme(legend.position="top")


a <- ggplot(AD_samplesheet,aes(BPC_UMAP_1,BPC_UMAP_2, col = celltype)) + geom_point(size = 0.5) + theme_minimal() + labs(x = "BinaryPC UMAP1",y = "BinaryPC UMAP2") +
   theme(legend.position="none")
b <- ggplot(AD_samplesheet,aes(PC_UMAP_1,PC_UMAP_2, col = celltype)) + geom_point(size = 0.5) + theme_minimal() + labs(x = "UMAP1",y = "UMAP2") +
  theme(legend.position="none")

a <- ggplot(AD_samplesheet,aes(PC_UMAP_1,PC_UMAP_2, col = as.factor(SYNPR_binary))) + geom_point(size = 0.25) + theme_minimal() +
  labs(x = "UMAP1",y = "UMAP2") + viridis::scale_color_viridis(discrete = T) + theme(legend.position="top")


b <- ggplot(AD_samplesheet,aes(PC_UMAP_1,PC_UMAP_2, col = SYNPR_normalized)) + geom_point(size = 0.25) + theme_minimal() +
  labs(x = "UMAP1",y = "UMAP2") + viridis::scale_color_viridis() + theme(legend.position="top")

a | b





