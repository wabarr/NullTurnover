#' Bin FADs from a tree, and test for pulses using chi-square
#' 
#' @param nBins the desired number of bins for combining FADs
#' @param myTree an object of class "phylo"
#' @param showBinPlot Whether or not to show the plots of binned residuals. Default is FALSE
#' @param showExpected whether or not to show the plot of the expected versus actual.  Default is FALSE
#' @param showTree Whether or not to plot the tree. Default is FALSE
#' 

DetectFADPulses <- function(nBins, myTree, showBinPlot = FALSE, showExpected = FALSE, showTree = FALSE){
  nls_mod <- FitExponential(myTree, showExpected = showExpected, showTree = showTree)
  MYA <- max(as.numeric(dist.nodes(myTree)))/2
  #set up the bounds for binning
  bounds <- seq(0, MYA, length.out = nBins + 1)
  lowerBounds <- bounds[1:length(bounds) - 1]
  upperBounds <- bounds[2:length(bounds)]
  
  sortedFADs <- GetFADs(myTree)
  expectedAddedTaxa <- predict(nls_mod, newdata=list(dates = upperBounds)) - predict(nls_mod, newdata=list(dates = lowerBounds))
  
  observedAddedTaxa <- sapply(upperBounds, FUN = function(x) sum(sortedFADs <= x)) - sapply(lowerBounds, FUN = function(x) sum(sortedFADs <= x))
  binnedResiduals <- observedAddedTaxa - expectedAddedTaxa
  
  mids <- (lowerBounds + upperBounds)/2
  if (showBinPlot){
    plot(mids, binnedResiduals, main = sprintf("%d taxa: %d bins (Observed - Expected)", length(myTree$tip.label), nBins))
    lines(mids, binnedResiduals, lty=2)
  }
  
  if (min(binnedResiduals, na.rm=TRUE)<0) warning("some intervals had fewer observed originations than expected from exponential model.  Running chi-square values only on intervals with positive values.")
  testResults <- chisq.test(binnedResiduals[binnedResiduals>0])
  print(testResults)
  return(testResults)
}