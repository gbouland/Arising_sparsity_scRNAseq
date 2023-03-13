library(harmony)
library(Seurat)
library(lisi)
## Binary
binary_brain <- readRDS("./Datasets_for_integration/To_big_for_github/brain_seurat_binary.rds")
binary_brain <- RunUMAP(binary_brain, dims = 1:10)

##Binary not integrated
DimPlot(binary_brain, reduction = "umap",group.by = "celltype")
DimPlot(binary_brain, reduction = "umap",group.by = "orig.ident")

##Binary integrated
binary_brain <- RunHarmony(binary_brain, "orig.ident")
binary_brain <- RunUMAP(binary_brain, reduction = "harmony", dims = 1:10)
DimPlot(binary_brain, reduction = "umap",group.by = "orig.ident")
DimPlot(binary_brain, reduction = "umap",group.by = "celltype")

##Binary LISI
LISI_binary <- compute_lisi(binary_brain@reductions$harmony@cell.embeddings[,1:10],binary_brain@meta.data,label_colnames = "orig.ident")
mean(LISI_binary$orig.ident)

## Counts
count_brain <- readRDS("./Datasets_for_integration/To_big_for_github/brain_seurat_normalized.rds")
count_brain <- RunUMAP(count_brain, dims = 1:10)

##count not integrated
DimPlot(count_brain, reduction = "umap",group.by = "celltype")
DimPlot(count_brain, reduction = "umap",group.by = "orig.ident")

##count integrated
count_brain <- RunHarmony(count_brain, "orig.ident")
count_brain <- RunUMAP(count_brain, reduction = "harmony", dims = 1:10)
DimPlot(count_brain, reduction = "umap",group.by = "orig.ident")
DimPlot(count_brain, reduction = "umap",group.by = "celltype")

##count LISI
LISI_count <- compute_lisi(count_brain@reductions$harmony@cell.embeddings[,1:10],count_brain@meta.data,label_colnames = "orig.ident")
mean(LISI_count$orig.ident)



