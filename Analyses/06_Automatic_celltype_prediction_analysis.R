files <- list.files("./To_big_for_github/Datasets_for_automatic_celltype_prediction/")
files <- data.frame(files = files, size = sapply(files,function(x)file.size(sprintf("./To_big_for_github/Datasets_for_automatic_celltype_prediction/%s",x))) |> unname())

runSingleR <- function(data, samplesheet, targetIDs, referencIDs){
  require(SingleR)
  predictions <- SingleR(test = data[,targetIDs],
                         ref = data[,referencIDs],
                         labels = samplesheet[match(referencIDs,samplesheet$colnames),"celltype"])
  return(as.factor(predictions$labels))
}

runScPred <- function(data, samplesheet, targetIDs, referencIDs){
  require(scPred)
  require(Seurat)
  reference <- data[,referencIDs]
  target <- data[,targetIDs]
  reference <- CreateSeuratObject(counts = reference)
  reference <- ScaleData(reference)
  reference <- RunPCA(reference,features = rownames(reference@assays$RNA@counts))
  target <- CreateSeuratObject(counts = target)
  target <- ScaleData(target)
  target <- RunPCA(target,features = rownames(reference@assays$RNA@counts))
  reference$celltype <- samplesheet[match(referencIDs,samplesheet$colnames),"celltype"]
  reference <- getFeatureSpace(reference, "celltype")
  reference <- trainModel(reference)
  target <- scPredict(target, reference)
  return(as.factor(target$scpred_no_rejection))
}

evaluate <- function(predicted, true){
  require(caret)
  to_eval <- data.frame("predicted" = predicted,
                        "true" = true)
  all_factors <- unique(c(levels(to_eval$predicted),levels(to_eval$true)))
  to_eval$predicted <- factor(to_eval$predicted,levels = all_factors)
  to_eval$true <- factor(to_eval$true,levels = all_factors)
  stats <- confusionMatrix(to_eval$predicted,to_eval$true)
  medianF1 <- median(stats$byClass[,7],na.rm = T)
  return(medianF1)
}

file <- "prep_AD_set.rds" #one of the 22 datasets
repl <- 1 #1-10

prepared_dataset <- readRDS(sprintf("./To_big_for_github/Datasets_for_automatic_celltype_prediction/%s",file))
samplesheet <- prepared_dataset[["samplesheet"]]
samplesheet$colnames <- as.character(samplesheet$colnames)
## Define reference and target datasets
referencIDs <- lapply(unique(samplesheet$celltype),function(celltype){
  celltype_colnames <- samplesheet[samplesheet$celltype == celltype,"colnames"]
  return(sample(celltype_colnames,size = round(length(celltype_colnames)*0.75)))
}) |> unlist()  
targetIDs <- samplesheet$colnames[!samplesheet$colnames %in% referencIDs]
true_labels <- as.factor(samplesheet[match(targetIDs,samplesheet$colnames),"celltype"])  
res <- lapply(c("runSingleR","runScPred"),function(runFun){
  message(runFun)
  runFun <- get(runFun)
  F1s <- lapply(c("normalized","binary","shuffled"),function(type){
    message(type)
    runFun(prepared_dataset[[type]], samplesheet, targetIDs, referencIDs)
  }) |> sapply(FUN = function(x){evaluate(predicted = x, true = true_labels)})
  names(F1s) <- c("normalized","binary","shuffled")
  return(F1s)
}) |> do.call(what = "rbind")
rownames(res) <- c("runSingleR","runScPred")
res <- as.data.frame(res)
filename <- file |>
  gsub(pattern = "prep_", replacement = "") |>
  gsub(pattern = ".rds", replacement = "")
res$dataset <- filename
res$no_of_celltypes <-  length(unique(true_labels))
res$smallest_celltype <- min(table(true_labels))
res$largest_celltype <- max(table(true_labels))


## Load Results and make plot
all_res <- lapply(list.files("./Data/automatic_celltype_prediction_results/"),
                  function(x)readRDS(sprintf("./Data/automatic_celltype_prediction_results/%s",x))) |> do.call(what = "rbind")
all_res$method <- rownames(all_res)
all_res$method <- gsub("([0-9]+).*$","",all_res$method)
all_res$method <- gsub("run","",all_res$method)
molten <- melt(all_res[,c("normalized","binary","shuffled","dataset","method")],id.vars = c("dataset","method"))
toMedian <- expand.grid(unique(molten$dataset),unique(molten$method),unique(molten$variable))
moltenMedian <- apply(toMedian,1,function(row){
  dataset <- row[1]
  method <- row[2]
  data_rep <-row[3]
  medianF1 <- median(molten[molten$dataset == dataset &
                              molten$method == method &
                              molten$variable == data_rep,"value"])
  variance <- var(molten[molten$dataset == dataset &
                           molten$method == method &
                           molten$variable == data_rep,"value"])
  data.frame("dataset" = dataset,
             "method" = method,
             "data_rep" = data_rep,
             "medianF1" = medianF1,
             "variance" = variance)
}) |> do.call(what = "rbind")


ggplot(moltenMedian,aes(x = dataset,y = data_rep, fill = medianF1)) + geom_tile() + viridis::scale_fill_viridis(limits = c(0,1)) +
  geom_text(aes(x = dataset,y = data_rep, label = round(medianF1,2)), color = "black", size = 4) + theme_minimal() + 
  facet_wrap(~method,nrow = 2) + scale_x_discrete(guide = guide_axis(angle = 45))
