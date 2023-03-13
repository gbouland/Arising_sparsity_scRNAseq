library(muscat)
library(scRNAseq)
library(scuttle)
library(magrittr)
library(ggplot2)
library(reshape2)

source("./Analyses/00_utils.R")



data(example_sce)
sce <- example_sce
ref <- prepSim(sce, verbose = FALSE)
pct_differentially_expressed <- c(0.5,0.4,0.3,0.2,0.1,0.01)
total_number_cells <- c(1000,5000,10000,50000)
number_of_samples_per_group <- c(5,10,20,50)
toRun_settings <- expand.grid(pct_differentially_expressed,
                              total_number_cells,
                              number_of_samples_per_group)

test <- lapply(1:nrow(toRun_settings),function(i){
  message(sprintf("main %s",i))
  pc_de <- toRun_settings[i,"Var1"]
  no_cells <- toRun_settings[i,"Var2"]
  no_samples <- toRun_settings[i,"Var3"]
  simulations <- lapply(1:10,function(x){
    message(x)
    simData(ref, p_dd = c(1-pc_de,0,pc_de,0,0,0),
            nc = no_cells, ng = 500, force = TRUE,nk = 1, ns = no_samples)
  })
  names(simulations) <- 1:length(simulations)
  res <- lapply(names(simulations),function(i){
    message(i)
    sim <- simulations[[i]]
    zeroAggr <- aggregate_data(sim,type = "zeros")
    meansAggr <- aggregate_data(sim,type = "mean")
    zerosRes <- runTests(zeroAggr,"ttest")
    meansRes <- runTests(meansAggr,"limmatrend")
    detectionRate <- rowMeans(sim@assays@data$counts>=1)
    true <- meansRes$true
    out <- data.frame("zerosDetected" = zerosRes$detected,
                      "meansDetected" = meansRes$detected,
                      "true" = zerosRes$true,
                      "detectionRate" = detectionRate,
                      "simRun" = i)
    return(out)
  }) |> do.call(what = rbind)
  res$pc_de <- pc_de
  res$no_cells <- no_cells
  res$no_samples <- no_samples
  return(res)
}) |> do.call(what = rbind)


test <- readRDS("./Data/comprehensive_simulation_results/comprehensive_simulation.rds")
test$groups <- cut(test$detectionRate, breaks = seq(from = 0, to = 1, by= 0.10))
test <- na.omit(test)
test$RUNID <- sprintf("%s:%s:%s:%s",test$simRun,test$pc_de,test$no_cells,test$no_samples)
uniqueIDs <- unique(test$RUNID)
plotData <- lapply(uniqueIDs,function(ID){
  zerosPred <- test[test$RUNID == ID,"zerosDetected"]
  meansPred <- test[test$RUNID == ID,"meansDetected"]
  true <- test[test$RUNID == ID,"true"]
  zerosF1 <- getMetric(y = true,x = zerosPred,metric = "FN")#metrics: TPR,FPR,PPV,FDR,FN,FP,F1
  meansF1 <- getMetric(y = true,x = meansPred,metric = "FN")#metrics: TPR,FPR,PPV,FDR,FN,FP,F1
  data.frame("zerosF1" = zerosF1,
             "meansF1" = meansF1,
             "simRun" = unlist(strsplit(ID,split = ":"))[1],
             "pc_de" = unlist(strsplit(ID,split = ":"))[2],
             "no_cells" = unlist(strsplit(ID,split = ":"))[3],
             "no_samples" = unlist(strsplit(ID,split = ":"))[4])
}) |> do.call(what = "rbind")

## boxPlots with simRun
plotData <- melt(plotData,id.vars = c("simRun","pc_de","no_cells","no_samples"))
plotData$no_cells <- factor(plotData$no_cells, levels = c("1000","5000","10000","50000"))
plotData$no_samples <- factor(plotData$no_samples, levels = c("5","10","20","50"))
ggplot(plotData,aes(x = pc_de ,value, col = variable)) + geom_boxplot() + facet_wrap(~no_samples+no_cells) + theme_minimal() + ylim(0,100)#+ ylim(0.5,1)




