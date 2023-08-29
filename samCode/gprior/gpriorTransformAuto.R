## g prior sampling using transform slice just one psuedo prior
## Author: Sam Johnson

source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# transform
gSamples <- gprior_sampler(list(Nsamples = 5000, Nburnin = 100,Nthin = 1, method = 'SteppingOut', w = 15))
# auto
autoList$approxSamples <- gSamples$g

# chainSamplesAuto <- sapply(1:nChains, FUN = \(i) { gprior_sampler(autoList) } )
samplesAuto <- gprior_sampler(autoList)

tempDf <- data.frame(nEval = sum(samplesAuto$nEval),
                     ESS = coda::effectiveSize(samplesAuto$g),
                     userTime = samplesAuto$time[1],
                     sysTime = samplesAuto$time[2],
                     elapsedTime = samplesAuto$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'auto',
                     mean = mean(samplesAuto$g),
                     lwr = quantile(samplesAuto$g, probs = 0.025),
                     upr = quantile(samplesAuto$g, probs = 0.975),
                     auc = auc(u = samplesAuto$u),
                     water = water_area(u = samplesAuto$u))
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/samplesAuto/samplesAuto',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')

