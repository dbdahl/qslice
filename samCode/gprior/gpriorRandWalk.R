## g prior sampling using rand walk
## Author: Sam Johnson

# setting up the sampler
source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# chainSamplesRandWalk <- sapply(1:nChains, FUN = \(i) { gprior_sampler(randWalkList) } )
samplesRandWalk <- gprior_sampler(randWalkList)

tempDf <- data.frame(nEval = mean(samplesRandWalk$nEval),
                     ESS = coda::effectiveSize(samplesRandWalk$g),
                     userTime = samplesRandWalk$time[1],
                     sysTime = samplesRandWalk$time[2],
                     elapsedTime = samplesRandWalk$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'randwalk',
                     mean = mean(samplesRandWalk$g),
                     lwr = quantile(samplesRandWalk$g, probs = 0.025),
                     upr = quantile(samplesRandWalk$g, probs = 0.975),
                     auc = 0,
                     water = 0)
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/randwalk/randwalk',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
