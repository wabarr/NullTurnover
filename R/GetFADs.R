#' Internal function for getting FADs from a tree
#' 
#' function to get ages (i.e. distance from root node) of each non-tip node
#' 
#' 
#' @param myTree an object of class "phylo"
#' 
GetFADs <- function(myTree){
  stopifnot(is(myTree, "phylo"))
  stopifnot(length(myTree$tip.label)>2)
  require(ape)
  
  allDistancesAsMatrix<-ape::dist.nodes(myTree)
  #tips are numbered 1 + nTips, so we can add 1 to get the root node
  rootNode<-length(myTree$tip.label)+1
  #we pull out the row for the root node of the distance matrix
  #looking only at columns for other nodes, not tips
  FADs<-allDistancesAsMatrix[rootNode,seq(from=rootNode+1,to=nrow(allDistancesAsMatrix))]   
  return(sort(as.vector(FADs)))
}