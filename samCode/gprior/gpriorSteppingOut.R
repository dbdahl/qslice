## g prior sampling using stepping out procedure
## Author: Sam Johnson

# setting up the sampler
source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# chainSamplesSteppingOut <- sapply(1:nChains, FUN = \(i) { gprior_sampler(steppingOutList) } )
samplesSteppingOut <- gprior_sampler(steppingOutList)

tempDf <- data.frame(nEval = sum(samplesSteppingOut$nEval),
                     ESS = coda::effectiveSize(samplesSteppingOut$g),
                     userTime = samplesSteppingOut$time[1],
                     sysTime = samplesSteppingOut$time[2],
                     elapsedTime = samplesSteppingOut$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'steppingout',
                     mean = mean(samplesSteppingOut$g),
                     lwr = quantile(samplesSteppingOut$g, probs = 0.025),
                     upr = quantile(samplesSteppingOut$g, probs = 0.975),
                     auc = 0,
                     water = 0)
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/steppingout/steppingout',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
