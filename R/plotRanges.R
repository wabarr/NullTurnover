#'
#'function plots fossil ranges using ggplot2
#
#'@param fossilRange A fossil range object produced by paleotree::fossilRecord2fossilRanges()
#'@param returnValue Do you want to simply return the plot, rather than printing it?

plotRanges <- function(fossilRange, returnValue=FALSE) {
  suppressPackageStartupMessages(require(ggplot2))
  theme_set(theme_bw(20))
  DF <- data.frame(fossilRange[(fossilRange[,1] > 0),]) #remove those with FAD of 0
  DF <- DF[complete.cases(DF),]
  DF <- DF[order(DF$FAD,decreasing = TRUE),]
  DF$num <- 1:nrow(DF)
  require(ggplot2)
  thePlot <- ggplot(DF, aes(x=FAD,xend=LAD, y=num, yend=num)) + 
    geom_segment(size=2, color="#045480") + 
    scale_x_reverse() + 
    labs(y="N (Cumulative Number of Taxa)", x="Time (Ma)")
  ifelse(returnValue,
         return(thePlot),
         print(thePlot))
}