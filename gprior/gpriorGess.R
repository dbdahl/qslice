## g prior sampling using gess slice sampler
## Author: Sam Johnson

# setting up the sampler
source('gpriorSetup.R')


# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      temp <- cucumber::slice_sampler_generalized_elliptical(x = g, log_target,
                                                             mu = 20, sigma = 15, df = 3, log = TRUE)
      g <- temp$x
      nEval <- temp$nEvaluations
      saving_updates()
    }
    burnin_thinning()
  })   
  save_time()
}


# saving the samples
saveRDS(chainSamples, file = 'data/gess.rds')

print('finished')