#' Test for pulses using using Foote Rates
#' 
#' @param treeDepth The temporal depth (in millions of years) represented by the tree.
#' @param criterion How extreme an observed turnover rate must be to be considered a turnover pulse (relative to the interquartile distance for the distribution of all bins).
#' @param nTaxa The desired number of taxa in the simulated phylogenetic tree. Passed to phytools::pbtree.
#' @param deathRate The rate of extinction.  Passed to phytools::pbtree.
#' @param desiredBinNumber The number of bins that treeDepth will be broken into to compute turnover rates. 
#' @param plotTree Whether or not to plot the simulated tree. Default is FALSE.
#' @param plotRates Whether or not to plot a histogram of the calculated rates.  Default is TRUE.
#' 
#' 
detectPulses <- function(treeDepth = 7, criterion = 1.5, nTaxa = 100, deathRate = 0, desiredBinNumber = 30, plotTree=FALSE, plotRates = FALSE) {
  require(phytools)
  require(paleotree)
  
  #TODO check on sample size varaiation
  myTree <- phytools::pbtree(n=nTaxa, d=deathRate, scale = treeDepth)
  
  ranges <- abs(phytools::nodeHeights(tree = myTree) - treeDepth)
  colnames(ranges) <- c("FAD", "LAD")
  
  binnedRanges <- paleotree::binTimeData(ranges, int.length = treeDepth / desiredBinNumber)
  
  if(plotTree) {
    plot(myTree) 
  }
  
  perCapRates <- as.data.frame(paleotree::perCapitaRates(binnedRanges, plot=FALSE))
  
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
  return(sum(perCapRates$pRate > quantile(perCapRates$pRate, 0.75, na.rm=T) * criterion, na.rm=TRUE))
  
}
