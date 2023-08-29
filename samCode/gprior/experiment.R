## expermiment



source('gpriorSetup.R')

library(ggplot2)


for(i in 1:5) {
  print(i)
  
  print('starting stepping out')
  # stepping out
  steppingOutList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
       method = 'SteppingOut', w = 25)
  
  steppingOutSamples <- replicate(1, gprior_sampler(steppingOutList))
  
  print('starting gess')
  # gess
  gessList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
       method = 'GESS', mu = 30, sigma = 50, df = 3)
  
  gessSamples <- replicate(1, gprior_sampler(gessList))
  
  print('starting latent')
  # latent
  latentList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
       method = 'Latent', s = 5, rate = 0.005)
  
  latentSamples <-replicate(1, gprior_sampler(latentList))
  
  print('starting randwalk')
  # randwalk
  randWalkList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
       method = 'RandWalk', c = 30)
  
  randWalkSamples <- replicate(1, gprior_sampler(randWalkList))
  
  print('starting indpendence')
  # indpendence
  burninSamples <- gprior_sampler(list(Nsamples = 5000, Nburnin = 100, Nthin = 1,
                                       method = 'SteppingOut', w = 25))
  
  independenceList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
                       method = 'Independence', pseudoType = 'OptimSamples',
                       approxSamples = burninSamples$g, everyIter = FALSE)
  
  independenceSamples <- replicate(1, gprior_sampler(independenceList))
  #replicate(1, gprior_sampler(independenceList))
  
  print('starting lapalce')
  # transform laplace
  laplaceList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
                    method = 'Transform', pseudoType = 'Laplace',
                    everyIter = 1)
  
  laplaceSamples <- replicate(1, gprior_sampler(laplaceList))
  
  print('starting auto')
  # tranform auto
  autoList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
                   method = 'Transform', pseudoType = 'Auto',
                   approxSamples = burninSamples$g, everyIter = FALSE)
  
  autoSamples <- replicate(1, gprior_sampler(autoList))
  
  print('starting optim samples')
  # transform optim samples
  optimSamplesList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
                           method = 'Transform', pseudoType = 'OptimSamples',
                           approxSamples = burninSamples$g, everyIter = FALSE)

  optimSampsSamples <- replicate(1, gprior_sampler(optimSamplesList))
  
  print('starting optim')

  # transform optim
  optimList <- list(Nsamples = 10000, Nburnin = 1000, Nthin = 2,
                    method = 'Transform', pseudoType = 'Optim', optimType = 'function',
                    everyIter = FALSE)

  optimSamples <- replicate(1, gprior_sampler(optimList))
  
  # combining the samples
  samplesList <- list(step = steppingOutSamples,
                      gess = gessSamples,
                      latent = latentSamples,
                      randWalk = randWalkSamples,
                      independence = independenceSamples,
                      lapalce = laplaceSamples,
                      auto = autoSamples,
                      optimSamps = optimSampsSamples,
                      optim = optimSamples)

  results <- lapply(samplesList, FUN = resultsFunc) |> plyr::ldply('data.frame')
  
  results |> 
    ggplot(aes(x = Est, y = .id, color = .id)) +
    geom_point() + 
    geom_errorbar(aes(xmin = lwrConf, xmax = uprConf)) +
    theme_minimal() +
    theme(
      legend.position = 'none'
    )
  
  ggsave(paste0('explorationPlots/MonteCarloErrorPlot',i,'.pdf'))
}

