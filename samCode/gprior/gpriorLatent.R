## g prior sampling using latent slice sampler
## Author: Sam Johnson

# setting up the sampler
source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# chainSamplesLatent <- sapply(1:nChains, FUN = \(i) { gprior_sampler(latentList) } )
samplesLatent <- gprior_sampler(latentList)

tempDf <- data.frame(nEval = sum(samplesLatent$nEval),
                     ESS = coda::effectiveSize(samplesLatent$g),
                     userTime = samplesLatent$time[1],
                     sysTime = samplesLatent$time[2],
                     elapsedTime = samplesLatent$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'latent',
                     mean = mean(samplesLatent$g),
                     lwr = quantile(samplesLatent$g, probs = 0.025),
                     upr = quantile(samplesLatent$g, probs = 0.975),
                     auc = 0,
                     water = 0)
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/latent/latent',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
