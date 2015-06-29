#' Internal Function to Fit exponential model to FADs
#' 
#' @param myTree an object of class "phylo"
#' @param showPlot whether or not to show plot.  Default is FALSE
#' 
#' @return a Nonlinear Least Squares object of class "nls" representing the exponential fit 

FitExponential <- function(myTree, showPlot = FALSE) {
  stopifnot(is(myTree, "phylo"))
  
  # assumes FADs are returned sorted from GetFADs() function
  sortedFADs <- GetFADs(myTree)
  
  N <- 1:length(sortedFADs)
  nls_mod <- nls(N ~ exp(a * sortedFADs), data=data.frame(sortedFADs = sortedFADs, N=N), start=list(a=1))
  
  if (showPlot){
    MYA <- max(as.numeric(dist.nodes(myTree)))/2
    #hypothetical date values for producing smooth exponential predictions
    evenFADs <- seq(0, MYA, length.out = 100)
    predictedN <- predict(nls_mod, newdata = sortedFADs)
    evenPredictedN <- predict(nls_mod, newdata=list(sortedFADs = evenFADs))
    plot(sortedFADs, N, pch=16,xlim=c(0,MYA))
    lines(evenFADs, evenPredictedN, col="red", lty=2, cex=10)
  }
  
  return(nls_mod)
}