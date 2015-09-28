#' Test for pulses using residuals from loess
#' 
#' @param nBins the desired number of bins for combining FADs
#' @param myTree an object of class "phylo". You must include EITHER a value for myTree OR a value for FADs, but not both.
#' @param FADs Vector of taxon first appearances. You must include EITHER a value for myTree OR a value for FADs, but not both.
#' @param showBinPlot Whether or not to show the plots of binned residuals. Default is FALSE
#' @param showExpected whether or not to show the plot of the expected versus actual.  Default is FALSE
#' @param showTree Whether or not to plot the tree. Default is FALSE
#' 

loessDetectPulses <- function(nBins, myTree=NULL, FADs=NULL, showPlots = FALSE, span = 0.2, showTree = FALSE){
  stopifnot(!(is.null(myTree) && is.null(FADs)))
  stopifnot(is.null(myTree) || is.null(FADs))
  require(phytools)
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
  bounds <- seq(round(minMYA,1), round(maxMYA,1), length.out = nBins + 1)
  lowerBounds <- bounds[1:length(bounds) - 1]
  upperBounds <- bounds[2:length(bounds)]
  mids <- (lowerBounds + upperBounds)/2
  
  FADindices <- 1:length(sortedFADs)
  loessMod <- loess(sortedFADs ~ FADindices, span = span)
  binGroups <- cut(sortedFADs, breaks=bounds)
  binnedResid <- tapply(resid(loessMod), INDEX=binGroups, FUN=sum)

  
  if (showPlots){
    par(mfrow=c(2,1))
    plot(y=FADindices, x=sortedFADs, pch=16, main
         =sprintf("Loess Line Fitted to FADs: span = %.2f", span))
    lines(y=FADindices, x=predict(loessMod, newdata = FADindices), col="red", lty=2, lwd=2)
    plot(x=mids, y=binnedResid, pch=16, main="Loess Residual Turnover Events in each bin")
    
  }
  
  return(
    data.frame(summedResid=scale(binnedResid), mids = mids)
    )
}