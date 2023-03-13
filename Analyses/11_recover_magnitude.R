library(ggplot2)

source("./BinaryMethods/functions/recoverMagnitude.R")

AD_set <- readRDS("./Data/prep_AD_set.rds")
AD_samplesheet <- AD_set$samplesheet
normalized <- AD_set$normalized

recovered <- recoverMagnitude(AD_set$binary)
plotData <- data.frame("recovered" = recovered["ADAP1",],
                       "normalized" = AD_set$normalized["ADAP1",])
plotData <- plotData[plotData$recovered!=0,]
ggplot(plotData,aes(recovered,normalized)) + geom_point()
