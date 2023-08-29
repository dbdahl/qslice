## g prior sampling using independence
## Author: Sam Johnson

# setting up the sampler
source('setupGprior.R')

i <- commandArgs(trailingOnly = TRUE)[1]

# samples to create the proposal disribution
gSamples <- gprior_sampler(list(Nsamples = 5000, Nburnin = 100,Nthin = 1, method = 'SteppingOut', w = 15))
# saving the samples 
independenceList$approxSamples <- gSamples$g

# chainSamplesIndependence <- sapply(1:nChains, FUN = \(i) { gprior_sampler(independenceList) } )
samplesIndependence <- gprior_sampler(independenceList)

tempDf <- data.frame(nEval = mean(samplesIndependence$nEval),
                     ESS = coda::effectiveSize(samplesIndependence$g),
                     userTime = samplesIndependence$time[1],
                     sysTime = samplesIndependence$time[2],
                     elapsedTime = samplesIndependence$time[3],
                     # sampPsec = metrics$EffSamp/metrics$userTime,
                     t = 'independence',
                     mean = mean(samplesIndependence$g),
                     lwr = quantile(samplesIndependence$g, probs = 0.025),
                     upr = quantile(samplesIndependence$g, probs = 0.975),
                     auc = 0,
                     water = 0)
tempDf$sampPsec <- tempDf$ESS/tempDf$userTime

# saving samples
write.table(tempDf, file = paste0('output/independence/independence',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
