#' Bin FADs from a tree, and test for pulses using chi-square, following Vrba (1995)
#' 
#' @param nBins the desired number of bins for combining FADs
#' @param myTree an object of class "phylo". You must include EITHER a value for myTree OR a value for FADs, but not both.
#' @param FADs Vector of taxon first appearances. You must include EITHER a value for myTree OR a value for FADs, but not both.
#' @param nls_mod Object of class 'nls' representing the exponential model fit between FADs and time.  If this is not provided, the model will be estimated from the FADs provided or the FADs inferred from the provided tree. 
#' @param showBinPlot Whether or not to show the plots of binned residuals. Default is FALSE
#' @param showExpected whether or not to show the plot of the expected versus actual.  Default is FALSE
#' @param showTree Whether or not to plot the tree. Default is FALSE
#' 

VrbaDetectPulses <- function(nBins, myTree=NULL, FADs=NULL, showBinPlot = FALSE, nls_mod = NULL, showExpected = FALSE, showTree = FALSE){
  stopifnot(!(is.null(myTree) && is.null(FADs)))
  stopifnot(is.null(myTree) || is.null(FADs))
  require(phytools)
  if (is.null(nls_mod)) nls_mod <- FitExponential(myTree = myTree, FADs = FADs, showExpected = showExpected, showTree = showTree)
  if (is.null(myTree)) {
    maxMYA <- max(FADs, na.rm=TRUE)
    minMYA <- min(FADs, na.rm=TRUE)
    sortedFADs <- sort(FADs)
  }
  else {
    maxMYA <- max(nodeHeights(myTree))
    minMYA <- 0
    sortedFADs <- GetFADs(myTree)
  }
  
  #set up the bounds for binning
  bounds <- seq(round(minMYA,1), round(maxMYA,1), by=0.1)
  lowerBounds <- bounds[1:length(bounds) - 1]
  upperBounds <- bounds[2:length(bounds)]
  
  expectedAddedTaxa <- predict(nls_mod, newdata=list(dates = upperBounds)) - predict(nls_mod, newdata=list(dates = lowerBounds))
  
  observedAddedTaxa <- sapply(upperBounds, FUN = function(x) sum(sortedFADs <= x)) - sapply(lowerBounds, FUN = function(x) sum(sortedFADs <= x))
  smallBinResiduals <- observedAddedTaxa - expectedAddedTaxa
  
  vrbaLargeBins <- gl(n=length(lowerBounds) / 3, k=3)
  
  binnedResiduals <- tapply(smallBinResiduals, INDEX=vrbaLargeBins, FUN=sum)
  mids <- (lowerBounds + upperBounds)/2
  binMids <- tapply(mids, INDEX=vrbaLargeBins, FUN=mean)
  
  if (showBinPlot){
    plot(binMids, binnedResiduals, main = sprintf("%d taxa: %d bins (Observed - Expected)\nExcluded Values in Red", length(sortedFADs), nBins))
    lines(binMids, binnedResiduals, lty=2)
    #points(binMids, binnedResiduals, col="red")
  }
  
  if (min(binnedResiduals, na.rm=TRUE)<0) warning("some intervals had fewer observed originations than expected from exponential model.  Running chi-square values only on intervals with positive values.")
  testResults <- chisq.test(binnedResiduals[binnedResiduals>0])
  print(testResults)
  return(list(chisq = testResults, lowerBounds = lowerBounds, upperBounds = upperBounds))
}