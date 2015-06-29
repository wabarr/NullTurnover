#' Internal function for getting FADs from a tree
#' 
#' function to get ages (i.e. distance from root node) of each non-tip node
#' 
#' 
#' @param myTree an object of class "phylo"
#' @param sorted
#' @param showTree Whether or not to plot the tree. Default is FALSE
#' 
GetFADs <- function(myTree, sorted = TRUE, showTree = FALSE){
  stopifnot(is(myTree, "phylo"))
  stopifnot(length(myTree$tip.label)>2)
  require(phytools)
    
  FADs <- phytools::nodeHeights(myTree)[,1]
  
  if (showTree) {
    plot(myTree,show.tip.label = FALSE, main=sprintf("Tree with %d tips", length(myTree$tip.label)))
  }
  
  ifelse(sorted, 
         return(sort(FADs)),
         return(FADs)
         )
}