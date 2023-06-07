## g prior sampling using latent slice sampler
## Author: Sam Johnson

# setting up the sampler
source('gpriorSetup.R')


# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      temp <- cucumber::slice_sampler_latent(x = g, s = 1, log_target,
                                             rate = 0.001, log = TRUE)
      g <- temp$x
      nEval <- temp$nEvaluations
      saving_updates()
    }
    burnin_thinning()
  })   
  save_time()
}


# saving samples
saveRDS(chainSamples, file = 'data/latent.rds')

print('finished')
