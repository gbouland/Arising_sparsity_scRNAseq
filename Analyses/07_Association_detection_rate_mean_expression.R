library(ggplot2)
detectionRates <- readRDS("./Data/detection_rate_and_mean_expression/MDD_PseudoDetectionRate.rds")
means <- readRDS("./Data/detection_rate_and_mean_expression/MDD_PseudoMean.rds")
#per individual
diag(cor(detectionRates,means,method = "spearman"))
plotData <- data.frame("Measurement rate" = detectionRates[,5],"Mean" = means[,5])
p1 <- ggplot(plotData,aes(y = Measurement.rate,x = log(Mean))) + geom_point() + theme_minimal()
p1
