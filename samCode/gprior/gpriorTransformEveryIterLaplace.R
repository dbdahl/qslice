## g prior sampling using transform slice psuedo updated at every interation using laplace
## Author: Sam Johnson

source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

laplaceList$everyIter <- 1

# chainSamplesLaplaceEveryIter <- sapply(1:nChains, FUN = \(i) { gprior_sampler(laplaceList) } )
samplesLaplaceEveryIter <- gprior_sampler(laplaceList)

tempDf <- data.frame(nEval = sum(samplesLaplaceEveryIter$nEval),
                     ESS = coda::effectiveSize(samplesLaplaceEveryIter$g),
                     userTime = samplesLaplaceEveryIter$time[1],
                     sysTime = samplesLaplaceEveryIter$time[2],
                     elapsedTime = samplesLaplaceEveryIter$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'laplace',
                     mean = mean(samplesLaplaceEveryIter$g),
                     lwr = quantile(samplesLaplaceEveryIter$g, probs = 0.025),
                     upr = quantile(samplesLaplaceEveryIter$g, probs = 0.975),
                     auc = auc(u = samplesLaplaceEveryIter$u),
                     water = water_area(u = samplesLaplaceEveryIter$u))
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/laplace/laplace',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')

