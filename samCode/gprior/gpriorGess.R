## g prior sampling using gess slice sampler
## Author: Sam Johnson

# setting up the sampler
source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# chainSamplesGESS <- sapply(1:nChains, FUN = \(i) { gprior_sampler(gessList) } )
samplesGESS <- gprior_sampler(gessList)

tempDf <- data.frame(nEval = sum(samplesGESS$nEval),
                     ESS = coda::effectiveSize(samplesGESS$g),
                     userTime = samplesGESS$time[1],
                     sysTime = samplesGESS$time[2],
                     elapsedTime = samplesGESS$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'gess',
                     mean = mean(samplesGESS$g),
                     lwr = quantile(samplesGESS$g, probs = 0.025),
                     upr = quantile(samplesGESS$g, probs = 0.975),
                     auc = 0,
                     water = 0)
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/gess/gess',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')

