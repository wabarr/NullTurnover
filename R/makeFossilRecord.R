#' Make fossil records
#' @param totalTime The temporal depth (in millions of years) represented by the simulation/
#' @param nTotalTaxa The minimum and maximum acceptable number of taxa in the tree.  A numeric vector of length two (e.g. c(min, max))
#' @param deathRate The average (per interval, per lineage) rate of extinction.  
#' @param birthRate The average (per interval, per lineage) rate of origination. 

makeFossilRecord <- function(birthRate, deathRate, totalTime = 7, nTotalTaxa=c(100,300)) {
  require(paleotree)
  fossRec <- simFossilRecord(p=birthRate, q=deathRate, totalTime = totalTime, nTotalTaxa = nTotalTaxa)
}
