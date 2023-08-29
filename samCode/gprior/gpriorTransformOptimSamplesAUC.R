## g prior sampling using transform slice just one psuedo prior found using optim based on samples
## Author: Sam Johnson

source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# transform
gSamples <- gprior_sampler(list(Nsamples = 5000, Nburnin = 100,Nthin = 1, method = 'SteppingOut', w = 15))
# optim samples
optimSamplesAUCList$approxSamples <- gSamples$g

# chainSamplesOptimSamples <- sapply(1:nChains, FUN = \(i) { gprior_sampler(optimSamplesAUCList) } )
samplesOptimSamplesAUC <- gprior_sampler(optimSamplesAUCList)

tempDf <- data.frame(nEval = sum(samplesOptimSamplesAUC$nEval),
                     ESS = coda::effectiveSize(samplesOptimSamplesAUC$g),
                     userTime = samplesOptimSamplesAUC$time[1],
                     sysTime = samplesOptimSamplesAUC$time[2],
                     elapsedTime = samplesOptimSamplesAUC$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'optimSamplesAUC',
                     mean = mean(samplesOptimSamplesAUC$g),
                     lwr = quantile(samplesOptimSamplesAUC$g, probs = 0.025),
                     upr = quantile(samplesOptimSamplesAUC$g, probs = 0.975),
                     auc = auc(u = samplesOptimSamplesAUC$u),
                     water = water_area(u = samplesOptimSamplesAUC$u))
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/optimSamplesAUC/optimSamplesAUC',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
