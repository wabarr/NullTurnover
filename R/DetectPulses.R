#' Test for turnover pulses in FADs
#' 
#' @param nBins the desired number of bins for combining FADs
#' @param myTree an object of class "phylo". You must include EITHER a value for myTree OR a value for FADs, but not both.
#' @param FADs Vector of taxon first appearances. You must include EITHER a value for myTree OR a value for FADs, but not both.
#' @param showBinPlot Whether or not to show the plots of binned residuals. Default is FALSE
#' @param showTree Whether or not to plot the tree. Default is FALSE
#' 
DetectPulses <- function(nBins, myTree=NULL, FADs=NULL, showBinPlot = FALSE, showTree = FALSE){
  
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
  mids <- (lowerBounds + upperBounds)/2
  
  observedAddedTaxa <- sapply(upperBounds, FUN = function(x) sum(sortedFADs <= x)) - sapply(lowerBounds, FUN = function(x) sum(sortedFADs <= x))
  plot(lowess(x = mids, y=observedAddedTaxa))
  }
