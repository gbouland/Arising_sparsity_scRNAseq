library(patchwork)
library(ggbeeswarm)
library(ggplot2)

files <- list.files("./Data/r_explained")
## Correlations between binarized representations and log-normalized
## representations of every cell, for every dataset
all_correlations <- lapply(files,function(file){
  res <- readRDS(sprintf("./Data/r_explained/%s",file)) |> as.data.frame()
  res$cor
}) |> unlist()
## Average over all cells
mean(min_cor,na.rm = T)
## Association of aforementioned correlation with detection rate (pct) and
## variance of non-zero counts 
out <- lapply(files,function(file){
  res <- readRDS(sprintf("./Data/r_explained/%s",file)) |> as.data.frame()
  res$comb <- res$pct * res$sd^2
  data.frame("pct" = cor(res$cor,res$pct,use = "complete.obs"),
             "sd" = cor(res$cor,res$sd,use = "complete.obs"),
             "pct_sd2" = cor(res$cor,res$comb,use = "complete.obs"))
}) |> do.call(what = "rbind")
## This was added to fullSummary.csv

## Figure 2a, and SFig. 1a,1b,1c (example: PaulHSCData.rds)
res <- readRDS(sprintf("./Data/r_explained/%s","PaulHSCData.rds")) |> 
  as.data.frame()
res$comb <- res$pct * res$sd^2

a <- ggplot(res,aes(cor,sd^2))  + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") +
  labs(x = "Correlation between binary and normalized cell",
       y = "Variance") +
  labs(title = "B")

b <- ggplot(res,aes(cor,pct))  + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") +
  labs(x = "Correlation between binary and normalized cell",
       y = "Detection rate") +
  labs(title = "A")

c <- ggplot(res,aes(cor,comb))  + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") +
  labs(x = "Correlation between binary and normalized cell",
       y = "Detection rate x Variance") +
  labs(title = "C")

b | a | c

## Figure 2b
summary_data <- read.csv2("./Data/fullSummary.csv")
abc <- lapply(unique(summary_data$Tech),function(x){
  median(summary_data[summary_data$Tech == x,"explainedByPCT_SD2"])
  }) |> unlist()
names(abc) <- unique(summary_data$Tech)
summary_data$Tech <- factor(summary_data$Tech,levels = names(sort(abc)))
d <- ggplot(summary_data[summary_data$explainedByPCT_SD2 < 0,],
            aes(Tech,explainedByPCT_SD2))  + 
  geom_boxplot()  + geom_quasirandom() +
  theme_minimal() + coord_flip()  + ylim(-1,0)

d




data <- readRDS(paste0("./Data/r_explained/",files[56]))
r_per_dataset <- lapply(files,function(x){
  data <- readRDS(paste0("./Data/r_explained/",x))
  data.frame(cor = data[,1],
             dataset = x)
}) |> do.call(what = "rbind")
r_per_dataset$dataset <- gsub(".rds","",r_per_dataset$dataset)
r_per_dataset$Tech <-  summary_data[match(r_per_dataset$dataset,summary_data$dataset),"Tech"]

library(ggridges)
ggplot(r_per_dataset,aes(x = cor,y = dataset, fill = Tech)) + geom_density_ridges() +
  xlim(0,1) + facet_grid(Tech ~ .,scales = "free_y", space = "free") + theme_minimal() + labs(x = "Correlation between binary and normalized expression (p)")







