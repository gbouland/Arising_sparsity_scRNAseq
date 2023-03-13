library(ggplot2)
summary_data <- read.csv2("./Data/fullSummary.csv")
##Correlation between year of publication and number of cells
cor.test(summary_data$year,summary_data$Number)
##Correlation between number of cells and detection rate
cor.test(summary_data$Number,summary_data$pct)
## Figure 1a
ggplot(summary_data, aes(year, Number)) +  scale_y_log10() +
  geom_smooth(method = "lm") + geom_point(aes(col = Tech)) + theme_minimal()
## Figure 2a
ggplot(summary_data, aes(Number,pct)) +  scale_x_log10() +
  geom_smooth(method = "lm") + geom_point(aes(col = Tech)) + theme_minimal()