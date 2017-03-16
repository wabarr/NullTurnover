# Analysis_MEE_manuscript_create_trees
Andrew Barr  
3/16/2017  


```r
#devtools::install_github("wabarr/NullTurnover")
library(NullTurnover)
library(paleotree)
```

```
## Loading required package: ape
```

```
## Warning: package 'ape' was built under R version 3.3.2
```



```r
iterations <- 1200
## Note....I am over producing fossil records because some of them will fail, and I need to get more than 1000
extinctionRates <- runif(iterations, 0.1, 0.6)
originationRates <- sapply(extinctionRates, FUN=function(x) runif(1, min = x, max=0.6))
treesOrErrors <- parallel::mclapply(1:iterations, mc.cores = 4, FUN=function(x){
  params <- list(
    origRateSim = originationRates[x],
    extRateSim = extinctionRates[x],
    totalTime=c(4,8),
    sampRateSim = runif(1, 0.5, 2),
    nruns = 1,
    plot=FALSE,
    startTaxa = round(runif(1, 1, 5), 0),
    nSamp = c(30,100),
    nExtant=c(0,500),
    nTotalTaxa = c(0,500),
    maxAttempts = 100
  )
  return(do.call(NullTurnover::simRecords,params)[[1]])
})

goodSimulations <- sapply(treesOrErrors, FUN=function(x) !is.na(x[1]))
goodSimulations <- treesOrErrors[goodSimulations]

saveRDS(sample(goodSimulations, 1000), "~/Dropbox/TurnoverPulseRedux/ThousandTrees.RDS")
```
