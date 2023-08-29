## g prior sampling using transform slice just one psuedo prior found using optim based on samples
## Author: Sam Johnson

source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# transform
gSamples <- gprior_sampler(list(Nsamples = 5000, Nburnin = 100,Nthin = 1, method = 'SteppingOut', w = 15))
# optim samples
optimSamplesList$approxSamples <- gSamples$g

samplesOptimSamples <- gprior_sampler(optimSamplesList)
# chainSamplesOptimSamples <- sapply(1:nChains, FUN = \(i) { gprior_sampler(optimSamplesList) } )

tempDf <- data.frame(nEval = sum(samplesOptimSamples$nEval),
                     ESS = coda::effectiveSize(samplesOptimSamples$g),
                     userTime = samplesOptimSamples$time[1],
                     sysTime = samplesOptimSamples$time[2],
                     elapsedTime = samplesOptimSamples$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'optimSamples',
                     mean = mean(samplesOptimSamples$g),
                     lwr = quantile(samplesOptimSamples$g, probs = 0.025),
                     upr = quantile(samplesOptimSamples$g, probs = 0.975),
                     auc = auc(u = samplesOptimSamples$u),
                     water = water_area(u = samplesOptimSamples$u))
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/optimSamples/optimSamples',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')

