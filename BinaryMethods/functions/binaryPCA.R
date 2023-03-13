##Normal PCA on binary single-cell data
##Using PCA function from Seurat
## binary   is binarized single cell matrix
## n        is number of components
binaryPCA <- function(binary, n = 10){
  require(Seurat)
  seurat <- CreateSeuratObject(counts = binary)
  seurat <- ScaleData(seurat, features = rownames(seurat))
  seurat <- RunPCA(seurat,features = rownames(seurat))
  binary_pca <- seurat@reductions$pca@cell.embeddings[,1:n]
  return(binary_pca)
}