library(BDA)
library(scRATE)
library(Seurat)

## M1 dataset ##
M1 <- readRDS("C:/Users/gabouland/Documents/004 PhD/zeros_review/data/M1.rds")
counts <- M1[[1]]
samplesheet <- M1[[2]]
rm(M1)
selection <- samplesheet[samplesheet$celltype %in% c("Exc L3-5 FEZF2 ASGR2", "Exc L3 THEMIS ENPEP"),]
table(selection$celltype)
counts <- counts[,selection$IDs]
object <- CreateSeuratObject(counts)
object[["celltype"]] <- selection$celltype
results <- BDA(object = object,
               contrast.by = "celltype", ident.1 = "Exc L3-5 FEZF2 ASGR2", ident.2 = "Exc L3 THEMIS ENPEP",
               p.adjust.method = "fdr", min.pct = 0.4, n.cores = 4, method = "LR")[[1]]

stable_genes <- rev(rownames(results[abs(results$Estimate) <= 0.025,]))[1:2]
top100 <- rownames(results[1:100,])

## Smart-seq dataset##
exons <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Allen_mtg/human_MTG_2018-06-14_exon-matrix.csv")
introns <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Allen_mtg/human_MTG_2018-06-14_intron-matrix.csv")
both <- exons + introns
rm(introns);rm(exons)
genes <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Allen_mtg/human_MTG_2018-06-14_genes-rows.csv")
rownames(both) <- genes$gene
both$V1 <- NULL
meta <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Allen_mtg/human_MTG_2018-06-14_samples-columns.csv")
meta <- meta[,c("sample_name","brain_subregion","class","donor","cluster")]
colnames(meta) <- c("IDs","celltype","status","patient","cluster")
selection <- meta[meta$patient == "H200.1030" & meta$cluster %in% c("Exc L4-5 RORB FOLH1B","Exc L4-6 RORB SEMA3E"),]
both <- both[,selection$IDs]

object <- CreateSeuratObject(both)
object[["celltype"]] <- selection$cluster
results <- BDA(object = object,
               contrast.by = "celltype", ident.1 = "Exc L4-5 RORB FOLH1B", ident.2 = "Exc L4-6 RORB SEMA3E",
               p.adjust.method = "fdr", min.pct = 0.5, n.cores = 4, method = "LR")[[1]]


stable_genes <- rev(rownames(results[abs(results$Estimate) <= 0.03,]))[1:100]
top100 <- rownames(results[1:100,])


stable_best_fit <- lapply(stable_genes,function(x){
  message(x)
  sel_counts <- as.numeric(counts[x,])
  gexpr <- data.frame(y=sel_counts, exposure=log(colSums(counts)))
  model_fit <- fit_count_models(gexpr,nCores = 1)
  elpd_loo <- compare_count_models(model_fit)
  best_fit <- rownames(elpd_loo)[1]
  message(best_fit)
  return(best_fit)
})
stable_fits <- data.frame("gene" = stable_genes, "best_fit" = unlist(stable_best_fit))

marker_best_fit <- lapply(top100,function(x){
  message(x)
  sel_counts <- as.numeric(counts[x,])
  gexpr <- data.frame(y=sel_counts, exposure=1)
  model_fit <- fit_count_models(gexpr,nCores = 1)
  elpd_loo <- compare_count_models(model_fit)
  best_fit <- rownames(elpd_loo)[1]
  message(best_fit)
  return(best_fit)
})
marker_fits <- data.frame("gene" = top100, "best_fit" = unlist(marker_best_fit))



## Fisher test ##
## 10x 
stable_10x <- readRDS("./Data/check_zero_inflation/topStable_res.rds")
difex_10x <- readRDS("./Data/check_zero_inflation/topDifEx_res.rds")
stable_10x$type <- "stable"
difex_10x$type <- "difex"
stable_10x$zi <- ifelse(grepl("ZI",stable_10x$res),"ZI","no_ZI")
difex_10x$zi <- ifelse(grepl("ZI",difex_10x$res),"ZI","no_ZI")
tenX <- rbind(stable_10x,difex_10x)
tenX$type <-as.factor(tenX$type)
tenX$zi <- as.factor(tenX$zi)
tenX$type <- factor(tenX$type,levels = c("stable","difex"))
tenXres <- fisher.test(table(tenX$type,tenX$zi))

## smart-seq
stable <- readRDS("./Data/check_zero_inflation/topStable_smart_seq.rds")
difex <- readRDS("./Data/check_zero_inflation/topDifEx_smart_seq.rds")
stable$type <- "stable"
difex$type <- "difex"
stable$zi <- ifelse(grepl("ZI",stable$res),"ZI","no_ZI")
difex$zi <- ifelse(grepl("ZI",difex$res),"ZI","no_ZI")
smartSeq <- rbind(stable,difex)
smartSeq$type <-as.factor(smartSeq$type)
smartSeq$zi <- as.factor(smartSeq$zi)
smartSeq$type <- factor(smartSeq$type,levels = c("stable","difex"))
smartSeqres <- fisher.test(table(smartSeq$type,smartSeq$zi))












