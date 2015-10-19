#' Test for pulses using using Foote Rates
#' 
#' @param fossilRecord Object of class fossilRecord, as produced by paleotree::makeFossilRecord()
#' @param totalTime The temporal depth (in millions of years) represented by the simulation/
#' @param criterion How extreme an observed turnover rate must be to be considered a turnover pulse (relative to the interquartile distance for the distribution of all bins).
#' @param desiredBinNumber The number of bins for computing turnover rates. \
#' @param targetBins Vector indicating subset of bins in which to look for pulses. You must provide only the interval starting time. 
#' @param plotTree Whether or not to plot the simulated tree. Default is FALSE.
#' @param plotRates Whether or not to plot a histogram of the calculated rates.  Default is TRUE.
#' 
#' 
detectPulses <- function(fossilRecord, criterion = 1.5, desiredBinNumber = 28, targetBins = seq(3.75, 1.5, by=-0.25), plotTree=FALSE, plotRates = FALSE) {
  require(paleotree)
  
  
  tryCatch({
    binnedRanges <- paleotree::binTimeData(fossilRecord2fossilTaxa(fossRec)[,c("orig.time", "ext.time")], int.length=totalTime/desiredBinNumber)
  },error = function(e){
    return(NA)
  })
    
  if(plotTree) {
    plot(myTree) 
  }
  
  perCapRates <- as.data.frame(paleotree::perCapitaRates(binnedRanges, plot=FALSE))
  perCapRates <- perCapRates[perCapRates[,1] %in% targetBins,]

  if(plotRates){
    require(ggplot2)
    thePlot <- qplot(x=(int.start + int.end)/2, y=pRate, data=perCapRates, geom="bar", stat="identity") + 
      geom_hline(yint=quantile(perCapRates$pRate, 0.75, na.rm=T) * criterion, color='red') + 
      scale_x_reverse() + 
      theme_bw(25) + 
      labs(x="Interval Midpoint (MYA)", y="Foote Origination Rate")
    print(thePlot)
  }
  
  # how many intervals fall outside the specified range compared to interquartile range
  pulses <- sum(perCapRates$pRate > quantile(perCapRates$pRate, 0.75, na.rm=T) * criterion, na.rm=TRUE)
  return(pulses)
  
}
