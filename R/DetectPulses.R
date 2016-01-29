#' Test for pulses using using Foote Rates
#' 
#' @param myTree simulated phylogenetic tree
#' @param criterion How extreme an observed turnover rate must be to be considered a turnover pulse (relative to the interquartile distance for the distribution of all bins).
#' @param desiredBinNumber The number of bins for computing turnover rates.
#' @param plotTree Whether or not to plot the simulated tree. Default is FALSE.
#' @param plotRates Whether or not to plot a histogram of the calculated rates.  Default is TRUE.
#' @param returnRatePlot Whether or not to add the ggplot of the rates to the return value.
#' 
detectPulses <- function(myTree, criterion = 1.5, desiredBinNumber = 10, plotTree=FALSE, plotRates = FALSE, returnRatePlot = FALSE) {
  require(phytools)
  require(paleotree)
  treeDepth <- max(diag(vcv.phylo(myTree)))
  ranges <- abs(phytools::nodeHeights(tree = myTree) - treeDepth)		
  colnames(ranges) <- c("FAD", "LAD")		
  
  binnedRanges <- paleotree::binTimeData(ranges, int.length = treeDepth / desiredBinNumber)
  if(plotTree) {
    plot(myTree) 
  }
  
  perCapRates <- as.data.frame(paleotree::perCapitaRates(binnedRanges, plot=FALSE))
  #get rid of first and last bin
  perCapRates <- perCapRates[-c(1, nrow(perCapRates)),]
  
  if(plotRates){
    require(ggplot2)
    thePlot <- 
      ggplot(data=perCapRates, aes(x=(int.start + int.end)/2, y=pRate)) + 
      geom_bar(stat="identity") + 
      geom_hline(yintercept = quantile(perCapRates$pRate, 0.75, na.rm=T) * criterion, color='red') + 
      scale_x_reverse() + 
      theme_bw(25) + 
      labs(x="Interval Midpoint (MYA)", y="Foote Origination Rate")
    print(thePlot)
  }
  
  # how many intervals fall outside the specified range compared to interquartile range
  pulses <- sum(perCapRates$pRate > quantile(perCapRates$pRate, 0.75, na.rm=T) * criterion, na.rm=TRUE)
  if(returnRatePlot) {
    return(list(pulses, thePlot))
  }
  return(pulses)
  
}






