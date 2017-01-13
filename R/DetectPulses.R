#' test for pulses using paleotree functions
#'

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
  # simulates multiple fossil records, outputs a list, the first element of which is theRecord
  # the second element of which is a dataframe showing the parameters used to simulate the record
    require(paleotree)
    # simulate taxon ranges under BDS model with simFossilRecord
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
          return(list(theRecord=NA, theParameters=theParameters))
          }
      )
    })
    return(fossilRecords)
}

detectPulses <- function(fossilRecords, binLength=0.25,criterion=1.5, plot=FALSE){
          ## takes an object of class fossilRecordPlusParams, produced by 
          ## simRecords() function
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



## plotRanges function plots fossil ranges using ggplot2
plotRanges <- function(fossilRange, returnValue=FALSE) {
  suppressPackageStartupMessages(require(ggplot2))
  theme_set(theme_bw(20))
  DF <- data.frame(fossilRange[(fossilRange[,1] > 0),]) #remove those with FAD of 0
  DF <- DF[complete.cases(DF),]
  DF <- DF[order(DF$FAD,decreasing = TRUE),]
  DF$num <- 1:nrow(DF)
  require(ggplot2)
  thePlot <- ggplot(DF, aes(x=FAD,xend=LAD, y=num, yend=num)) + 
    geom_segment(size=2, color="#045480") + 
    scale_x_reverse() + 
    labs(y="N (Cumulative Number of Taxa)", x="Time (Ma)")
  ifelse(returnValue,
         return(thePlot),
         print(thePlot))
}
