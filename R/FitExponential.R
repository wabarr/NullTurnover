#' Internal Function to Fit exponential model to FADs
#' Function is to recreate images from Vrba (1995) classic turnover pulse figure
#' 
#' @param myTree an object of class "phylo". You must include EITHER a value for myTree OR a value for FADs
#' @param showExpected whether or not to show the plot of the expected versus actual.  Default is FALSE
#' @param showTree Whether or not to plot the tree. Default is FALSE
#' @param FADs Vector of taxon first appearances.  You must include EITHER a value for myTree OR a value for FADs
#' @return a Nonlinear Least Squares object of class "nls" representing the exponential fit 

fitExponential <- function(myTree=NULL, FADs=NULL, showExpected = FALSE, showTree = FALSE) {
  stopifnot(!(is.null(myTree) && is.null(FADs)))
  stopifnot(is.null(myTree) || is.null(FADs))
  
  if (is.null(myTree)) sortedFADs <- sort(FADs)
    else sortedFADs <- GetFADs(myTree, showTree = showTree) # assumes FADs are returned sorted from GetFADs() function
  
  N <- 1:length(sortedFADs)
  nls_mod <- nls(N ~ exp(a * dates), data=data.frame(dates = sortedFADs, N=N), start=list(a=1))
  
  if (showExpected){
    MYA <- max(nodeHeights(myTree))
    #hypothetical date values for producing smooth exponential predictions
    evenFADs <- seq(0, MYA, length.out = 100)
    predictedN <- predict(nls_mod, newdata = sortedFADs)
    evenPredictedN <- predict(nls_mod, newdata=list(dates = evenFADs))
    plot(sortedFADs, N, pch=16,xlim=c(0,MYA), main=sprintf("%d taxa", length(myTree$tip.label)))
    lines(evenFADs, evenPredictedN, col="red", lty=2, cex=10)
  }
  
  return(nls_mod)
}