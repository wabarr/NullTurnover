#' Bin FADs from a tree, and test for pulses using chi-square
#' 
#' @param nBins the desired number of bins for combining FADs
#' @param myTree an object of class "phylo"
#' 
#' 

TestFADs <- function(nBins, myTree){
  nls_mod <- FitExponential(myTree)
  MYA <- max(as.numeric(dist.nodes(myTree)))/2
  #set up the bounds for binning
  bounds <- seq(0, MYA, length.out = nBins + 1)
  lowerBounds <- bounds[1:length(bounds) - 1]
  upperBounds <- bounds[2:length(bounds)]
  
  sortedFADs <- GetFADs(myTree)
  expectedAddedTaxa <- predict(nls_mod, newdata=list(sortedFADs = upperBounds)) - predict(nls_mod, newdata=list(sortedFADs = lowerBounds))
  
  observedAddedTaxa <- sapply(upperBounds, FUN = function(x) sum(sortedFADs <= x)) - sapply(lowerBounds, FUN = function(x) sum(sortedFADs <= x))
  binnedResiduals <- observedAddedTaxa - expectedAddedTaxa
  
  mids <- (lowerBounds + upperBounds)/2
  plot(mids, binnedResiduals, main = "Observed new FADs - Expected ")
  lines(mids, binnedResiduals, lty=2)
  
  return(chisq.test(binnedResiduals))
}