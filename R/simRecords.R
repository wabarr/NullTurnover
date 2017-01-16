#' simulate fossil records.  This function is a wrapper for paleotree::simFossilRecord that does error catching
#'@return a list of class fossilRecordPlusParams which includes the simulated record (or an error), plus a list of parameters used in the simulation.

simRecords <- function(origRateSim, 
                       extRateSim, 
                       nruns,
                       totalTime, 
                       sampRateSim,
                       startTaxa,
                       maxAttempts,
                       nSamp=c(0,500),
                       nExtant=c(10,500),
                       nTotalTaxa=c(10,500),
                       plot=TRUE,
                       print.runs=FALSE) {

  require(paleotree)
  # simulate taxon ranges under BD model with simFossilRecord
  fossilRecords<-lapply(1:nruns,FUN=function(x) {
    tryCatch({
      theRecord <- simFossilRecord(p=origRateSim,
                                   q=extRateSim, 
                                   r = sampRateSim,
                                   totalTime = totalTime,
                                   nTotalTaxa = nTotalTaxa,
                                   nExtant = nExtant,
                                   nSamp = nSamp,
                                   startTaxa = startTaxa,
                                   print.runs = print.runs,
                                   plot = plot,
                                   maxAttempts=maxAttempts
      )
      theParameters <- data.frame(origRateSim,extRateSim,sampRateSim,totalTime,nTotalTaxa,nExtant,nSamp,startTaxa)
      returnVal <- list(theRecord=theRecord,theParameters=theParameters)
      class(returnVal) <- c(class(returnVal), "fossilRecordPlusParams")
      return(returnVal)
    }
    ,error=function(e){
      theParameters <- data.frame(origRateSim,extRateSim,sampRateSim,totalTime,nTotalTaxa,nExtant,nSamp,startTaxa)
      returnVal <- list(theRecord=NA, theParameters=theParameters)
      class(returnVal) <- c(class(returnVal), "fossilRecordPlusParams")
      return(returnVal)
    }
    )
  })
  return(fossilRecords)
}