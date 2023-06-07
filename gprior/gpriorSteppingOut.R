## g prior sampling using stepping out procedure
## Author: Sam Johnson

# setting up the sampler
source('gpriorSetup.R')


# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      temp <- cucumber::slice_sampler_stepping_out(g, log_target,
                                                   w = 50, log = TRUE, max = Inf)
      g <- temp$x
      nEval <- temp$nEvaluations
      saving_updates()
    }
    burnin_thinning()
  })   
  save_time()
}

# saving samples
saveRDS(chainSamples, file = 'data/steppingOut.rds')

print('finished')
