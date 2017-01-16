#' test for pulses using paleotree functions
#' @param fossilRecords  an object of class fossilRecordPlusParams, produced by NullTurnover::simRecords() function
#' @param binLength Length of time corresponding to bins, in simulation time units. Default is 0.25
#' @param criterion Criterion for determining if a bin is a turnover pulse or not. Based on upper 3rd quartile of turnover rates for each bin.  Default 1.5 means that a pulse is defined as 1.5 times higher than the 3rd quartile, including zeros. 
#' @param plot Whether or not to plot the ranges using NullTurnover::plotRanges(). Defaults to FALSE
detectPulses <- function(fossilRecords, binLength=0.25,criterion=1.5, plot=FALSE){

          stopifnot(inherits(fossilRecords[[1]], "fossilRecordPlusParams"))
          require(paleotree)
          require(ggplot2)
          # get ranges
          fossilRanges<-lapply(fossilRecords,FUN = function(x) {fossilRecord2fossilRanges(x$theRecord)})
          # get orig rates
          origRateList<-lapply(fossilRanges,getOrigRates, binLength=binLength)
          
          # test for pulses
          pulseTestList<-lapply(origRateList,testPulses,criterion=criterion)
          # return 
          
          ifelse(plot,{
            invisible(lapply(1:length(fossilRanges), FUN=function(x){
              thePlot <- plotRanges(fossilRanges[[x]], returnValue=TRUE)
              print(pulseTestList[[x]])
              if(any(pulseTestList[[x]]$isPulse, na.rm=TRUE)){ #if there are pulses, add the boxes to the plot
                thePlot <- thePlot + geom_rect(data=subset(pulseTestList[[x]], isPulse==TRUE),
                                               inherit.aes=FALSE,
                                               #ymax has to be set outsid aes() because it isn't mapped to a column
                                               ymax=sum(fossilRanges[[x]][,1] > 0, na.rm = T), #how many FADs are > 0 and not NA
                                               mapping=aes(xmin=int.start,xmax=int.end,ymin=0),
                                               fill="transparent",
                                               color="red",
                                               linetype="dashed")
              }
              print(thePlot)
              return(pulseTestList)
            }))
          },
          return(pulseTestList))
          
}


getOrigRates<-function(ranges,binLength){
  ## process a taxon ranges to get footes origination rates
  # bin the data
  binnedRanges <- binTimeData(ranges, int.length = binLength)
  # get foote rates
  perCapRates <- perCapitaRates(binnedRanges,plot=FALSE)
  # get rid of first and last bin
  origRates <- perCapRates[-c(1, nrow(perCapRates)),c("int.start", "int.end", "pRate")]
  return(origRates)
}

testPulses<-function(origRate,criterion){
  #process the output of getOrigRates() to test for pulses
  origRate <- data.frame(origRate)
  # how many intervals fall outside the specified range
  # compared to interquartile range ? (Bibi and Kiessling 2015)
  IQRcrit<-quantile(origRate[,"pRate"], 0.75, na.rm=TRUE) * criterion
  origRate$isPulse <- origRate[,"pRate"] > IQRcrit
  return(origRate)
}


