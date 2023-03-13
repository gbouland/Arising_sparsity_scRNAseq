aggregate_zeros <- function(sceObject){
  counts <- sceObject@assays@data$counts
  gi <- metadata(sceObject)$gene_info
  cellsheet <- data.frame("colnames" = colnames(counts),
                          "group" = sceObject@colData@listData$group_id,
                          "sampleID" = as.character(sceObject@colData@listData$sample_id))
  rownames(counts) <- paste0(rownames(counts),"-",gi$category)
  samplesIDs <- unique(cellsheet$sampleID)
  zero_aggr <- lapply(samplesIDs,function(ID){
    cells <- cellsheet[cellsheet$sampleID == ID,"colnames"]
    rowMeans(data.frame(counts[,cells] >=1))
  }) %>% do.call(what = cbind)
  colnames(zero_aggr) <- samplesIDs
  return(zero_aggr)
}

aggregate_mean <- function(sceObject){
  res <- aggregateData(sceObject,by = "sample_id",fun = "mean")
  aggregated <- res@assays@data[[1]]
  rownames(aggregated) <- paste0(rownames(aggregated),"-",metadata(res)$gene_info$category)
  return(aggregated)
}

aggregate_sum <- function(sceObject){
  res <- aggregateData(sceObject,by = "sample_id",fun = "sum")
  aggregated <- res@assays@data[[1]]
  rownames(aggregated) <- paste0(rownames(aggregated),"-",metadata(res)$gene_info$category)
  return(aggregated)
}

aggregate_data <- function(sceObject, type){
  aggregated <- switch(type,
                       "zeros" = aggregate_zeros(sceObject),
                       "mean" = aggregate_mean(sceObject),
                       "sum" = aggregate_sum(sceObject))
  return(aggregated)
}


runTests <- function(Data, test){
  res <- switch(test,
                "limmavoom" = runLimma(Data),
                "limmatrend"= runLimmaTrend(Data),
                "edgeR" = runEdgeR(Data),
                "ttest" = runTtest(Data),
                "wilcox" = runWilcox(Data))
  
  
  res$true <- as.numeric(grepl("-de",rownames(res)))
  res$detected <- ifelse(res$fdr <= 0.05,1,0)
  return(res)
}

runLimmaTrend <- function(Data){
  require(limma)
  require(edgeR)
  groups <- strsplit(colnames(Data),split = "[.]") %>% sapply(FUN = function(x)x[2])
  dge <- DGEList(Data, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~groups)
  logCPM <- cpm(dge, log=FALSE, prior.count=0.1)
  fit <- lmFit(logCPM, design = design)
  fit <- eBayes(fit,trend = TRUE)
  res <- data.frame("pvalue" = fit$p.value[,2])
  res$fdr <- p.adjust(res$pvalue,method = "BH")
  return(res)
}

runLimma <- function(Data){
  require(limma)
  require(edgeR)
  groups <- strsplit(colnames(Data),split = "[.]") %>% sapply(FUN = function(x)x[2])
  dge <- DGEList(Data, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~groups)
  vm <- voom(dge, design = design)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  res <- data.frame("pvalue" = fit$p.value[,2])
  res$fdr <- p.adjust(res$pvalue,method = "BH")
  return(res)
}

runEdgeR <- function(Data){
  require(limma)
  require(edgeR)
  groups <- strsplit(colnames(Data),split = "[.]") %>% sapply(FUN = function(x)x[2])
  dge <- DGEList(Data, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~groups)
  dge <- estimateDisp(dge, design = design)
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  res <- topTags(lrt, n = nrow(Data))$table
  res$fdr <- p.adjust(res$PValue,method = "BH")
  return(res)
}


runTtest <- function(Data){
  res <- apply(Data,1,function(row){
    if(length(unique(row)) == 1){
      return(data.frame(t = NA, p = 1))
    }else{
      res_t <- t.test(row[grepl(".A",names(row))],row[grepl(".B",names(row))])
      data.frame(t = res_t$statistic, p = res_t$p.value)
    }
  }) %>% do.call(what = rbind)
  res$fdr <- p.adjust(res$p,method = "BH")
  return(res)
}


runWilcox <- function(Data){
  index <- grep(".A",colnames(Data))
  wilcox_P <- sapply(rownames(Data),function(gene){
    min(2 * min(limma::rankSumTestWithCorrelation(index = index, statistics = Data[gene,])), 1)
  })
  res <- data.frame("P" = wilcox_P, "fdr" = p.adjust(wilcox_P,method = "BH"))
  return(res)
  
}
getMetric <- function(y, x , metric){
  tp <- sum(y == x & y == 1)
  fp <- sum(x == 1 & y == 0)
  fn <- sum(x == 0 & y == 1)
  tn <- sum(x == y & y == 0)
  switch(metric,
         "TPR" = {
           tp / (tp +fn)
         },
         "FPR" = {
           fp / (fp + tn)
         },
         "PPV" = {
           tp / (tp +fp)
         },
         "FDR" = {
           fp / (tp +fp)
         },
         "FN" = {
           fn
         },
         "FP" = {
           fp
         },
         "F1" = {
           tp  / (tp + 0.5*(fp+fn))
         })
}