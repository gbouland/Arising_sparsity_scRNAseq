library(reshape2)
library(patchwork)
#library(scBFA)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(cluster)
library(locStra)
library(Matrix)
library(svd)
library(uwot)
AD_set <- readRDS("./Data/prep_AD_set.rds")
AD_samplesheet <- AD_set$samplesheet
## scBFA
sce <- SingleCellExperiment(assay = list(logcounts = AD_set$binary))
bfa_model = scBFA(scData = sce, numFactors = 10)
cell_loadings <- bfa_model$ZZ
rownames(cell_loadings) <- colnames(AD_set$normalized)
colnames(cell_loadings) <- paste0("scBFA_",1:10)
saveRDS(cell_loadings, "./Data/reducedDimensions/AD_set_scBFA_10.rds")
##Jaccard Eigen vectors
sparseM <- Matrix(AD_set$binary,sparse=TRUE)
jac <- jaccardMatrix(sparseM, useCpp = TRUE, sparse = TRUE)
jac_eigen <- trlan.eigen(jac,neig = 10)
jac_eigen <- jac_eigen$u
rownames(jac_eigen) <- colnames(AD_set$normalized)
colnames(jac_eigen) <- paste0("Jaccard_",1:10)
saveRDS(jac_eigen, "./Data/reducedDimensions/AD_set_JaccardEigen_10.rds")
##Binary PCA
AD_seurat <- CreateSeuratObject(counts = AD_set$binary, project = "AD")
AD_seurat <- ScaleData(AD_seurat, features = rownames(AD_seurat))
AD_seurat <- RunPCA(AD_seurat,features = rownames(AD_seurat))
binary_pca <- AD_seurat@reductions$pca@cell.embeddings[,1:10]
saveRDS(binary_pca, "./Data/reducedDimensions/AD_set_BinaryPCA_10.rds")
## Count PCA
AD_seurat <- CreateSeuratObject(counts = AD_set$normalized, project = "AD")
AD_seurat <- ScaleData(AD_seurat, features = rownames(AD_seurat))
AD_seurat <- RunPCA(AD_seurat,features = rownames(AD_seurat))
normal_pca <- AD_seurat@reductions$pca@cell.embeddings[,1:10]
#saveRDS(normal_pca, "./Data/reducedDimensions/AD_set_NormalPCA_10.rds")
data <- c("./Data/reducedDimensions/AD_set_scBFA_10.rds",
          "./Data/reducedDimensions/AD_set_JaccardEigen_10.rds",
          "./Data/reducedDimensions/AD_set_BinaryPCA_10.rds",
          "./Data/reducedDimensions/AD_set_NormalPCA_10.rds")
data <- lapply(data, readRDS)
names(data) <- c("scBFA","Jaccard","BinaryPCA","NormalPCA")
AD_samplesheet$BF1 <- data$scBFA[,1]
AD_samplesheet$BF2 <- data$scBFA[,2]
AD_samplesheet$JEV1 <- data$Jaccard[,1]
AD_samplesheet$JEV2 <- data$Jaccard[,2]
AD_samplesheet$BPC1 <- data$BinaryPCA[,1]
AD_samplesheet$BPC2 <- data$BinaryPCA[,2]
AD_samplesheet$PC1 <- data$NormalPCA[,1]
AD_samplesheet$PC2 <- data$NormalPCA[,2]
BF <- ggplot(AD_samplesheet,aes(BF1,BF2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() +
  theme(legend.position = "none") + labs(title = "A")
JEV <- ggplot(AD_samplesheet,aes(JEV1,JEV2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() +
  theme(legend.position = "none") + labs(title = "B")
BPC <- ggplot(AD_samplesheet,aes(BPC1,BPC2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() +
  theme(legend.position = "none") + labs(title = "C")
PC <- ggplot(AD_samplesheet,aes(PC1,PC2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() +
  theme(legend.position = "none") + labs(title = "D")
BF | JEV | BPC | PC
BF_UMAP <- umap(data$scBFA,n_neighbors = 10, min_dist= 0.3,spread = 1)
JEV_UMAP <- umap(data$Jaccard,n_neighbors = 10, min_dist= 0.3,spread = 1)
BPC_UMAP <- umap(data$BinaryPCA,n_neighbors = 10, min_dist= 0.3,spread = 1)
PC_UMAP <- umap(data$NormalPCA,n_neighbors = 10, min_dist= 0.3,spread = 1)
AD_samplesheet$BF_UMAP_1 <- BF_UMAP[,1]
AD_samplesheet$BF_UMAP_2 <- BF_UMAP[,2]
AD_samplesheet$JEV_UMAP_1 <- JEV_UMAP[,1]
AD_samplesheet$JEV_UMAP_2 <- JEV_UMAP[,2]
AD_samplesheet$BPC_UMAP_1 <- BPC_UMAP[,1]
AD_samplesheet$BPC_UMAP_2 <- BPC_UMAP[,2]
AD_samplesheet$PC_UMAP_1 <- PC_UMAP[,1]
AD_samplesheet$PC_UMAP_2 <- PC_UMAP[,2]
UMAP_BF <- ggplot(AD_samplesheet,aes(BF_UMAP_1,BF_UMAP_2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() + theme(legend.position = "none") +
  labs(title = "A")
UMAP_JEV <- ggplot(AD_samplesheet,aes(JEV_UMAP_1,JEV_UMAP_2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() + theme(legend.position = "none") +
  labs(title = "B")
UMAP_BPC <- ggplot(AD_samplesheet,aes(BPC_UMAP_1,BPC_UMAP_2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() + theme(legend.position = "none") +
  labs(title = "C")
UMAP_PC <- ggplot(AD_samplesheet,aes(PC_UMAP_1,PC_UMAP_2, col = celltype)) +
  geom_point(size = 0.15) + theme_minimal() + theme(legend.position = "none") +
  labs(title = "D")
UMAP_BF | UMAP_JEV | UMAP_BPC | UMAP_PC
## Silhoutte
AD_samplesheet$clusterNUM <- as.numeric(as.factor(AD_samplesheet$celltype))
sill <- lapply(data,function(dat){
  comps <- dat
  dists <- dist(comps,method = "euclidean")
  sill <- cluster::silhouette(x = AD_samplesheet$clusterNUM, dist = dists)
  summ <- summary(sill)
  out <- c(summ$clus.avg.widths,summ$avg.width)
  names(out) <- c(levels(AD_samplesheet$celltype),"average")
  return(out)
}) |> do.call(what = "cbind")
colnames(sill) <- c("scBFA","JEVs","BinaryPCA","countPCA")
sill_umap <- lapply(data,function(dat){
  comps <- dat
  UMAP <- uwot::umap(comps,n_neighbors = 10, min_dist= 0.3,spread = 1, 
                     n_components = 2)
  dists <- dist(UMAP,method = "euclidean")
  sill <- cluster::silhouette(x = AD_samplesheet$clusterNUM, dist = dists)
  summ <- summary(sill)
  out <- c(summ$clus.avg.widths,summ$avg.width)
  names(out) <- c(levels(AD_samplesheet$celltype),"average")
  return(out)
}) |> do.call(what = "cbind")
colnames(sill_umap) <- c("scBFA","JEVs","BinaryPCA","countPCA")
molten_sill <- melt(sill)
molten_sill_umap <- melt(sill_umap)
molten_sill$type <- "Reduced dimensions (n=10)"
molten_sill_umap$type <- "UMAP dimensions (n=2)"
molten <- rbind(molten_sill,molten_sill_umap)
ggplot(molten,aes(Var1,Var2,fill = value)) + geom_tile() + facet_grid(~type) +
  viridis::scale_fill_viridis() + theme_minimal() + 
  labs(x = "cell type",y = "method", fill ="Average silhouette\nwidth") +
  geom_text(aes(Var1, Var2, label = round(value,2)), color = "white", size = 4)




## pairwise distance
data <- c("./Data/reducedDimensions/AD_set_scBFA_10.rds",
          "./Data/reducedDimensions/AD_set_JaccardEigen_10.rds",
          "./Data/reducedDimensions/AD_set_BinaryPCA_10.rds",
          "./Data/reducedDimensions/AD_set_NormalPCA_10.rds")
data <- lapply(data, readRDS)
names(data) <- c("scBFA","Jaccard","BinaryPCA","NormalPCA")
BF_UMAP <- umap(data$scBFA,n_neighbors = 10, min_dist= 0.3,spread = 1)
JEV_UMAP <- umap(data$Jaccard,n_neighbors = 10, min_dist= 0.3,spread = 1)
BPC_UMAP <- umap(data$BinaryPCA,n_neighbors = 10, min_dist= 0.3,spread = 1)
PC_UMAP <- umap(data$NormalPCA,n_neighbors = 10, min_dist= 0.3,spread = 1)


sel <- sample(rownames(data$BinaryPCA),5000)
BPC_UMAP <- BPC_UMAP[sel,]
PC_UMAP <- PC_UMAP[sel,]
BF_UMAP <- BF_UMAP[sel,]
JEV_UMAP <- JEV_UMAP[sel,]

UMAPs <- list(BPC_UMAP,
              PC_UMAP,
              BF_UMAP,
              JEV_UMAP)
all_pairwise_dists <- lapply(UMAPs,function(UMAP){
  message("start")
  pairwise_dist <- as.matrix(dist(UMAP))
  pairwise_dist[upper.tri(pairwise_dist)] <- NA
  diag(pairwise_dist) <- NA
  pairwise_dist <- melt(pairwise_dist)
  pairwise_dist <- pairwise_dist[!is.na(pairwise_dist$value),]
  return(pairwise_dist$value)
}) |> do.call(what = "cbind")
colnames(all_pairwise_dists) <- c("binaryPCA","countPCA","scBFA","JaccardEV")
all_pairwise_dists <- as.data.frame(all_pairwise_dists)
plotData <- all_pairwise_dists[sample(1:nrow(all_pairwise_dists),10000),]
plotData <- as.data.frame(plotData)
plotData$bins <- cut(plotData$countPCA,breaks = seq(0,24,by = 1))
cor(all_pairwise_dists$countPCA,all_pairwise_dists$binaryPCA,method = "pearson")
cor(all_pairwise_dists$countPCA,all_pairwise_dists$scBFA,method = "pearson")
cor(all_pairwise_dists$countPCA,all_pairwise_dists$JaccardEV,method = "pearson")
a <- ggplot(plotData,aes(bins,binaryPCA)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "binary-PCA based UMAP pairwise distances", x = "binned count-based UMAP pairwise distances", title = "A") + theme(plot.title = element_text(size=22))

b <- ggplot(plotData,aes(bins,JaccardEV)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "JaccardEV based UMAP pairwise distances", x = "binned count-based UMAP pairwise distances", title = "B") + theme(plot.title = element_text(size=22))

c <- ggplot(plotData,aes(bins,scBFA)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "scBFA based UMAP pairwise distances", x = "binned count-based UMAP pairwise distances", title = "C") + theme(plot.title = element_text(size=22))


a | b | c 


