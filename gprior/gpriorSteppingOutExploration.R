## exploring stepping out for the g prior
## Author: Sam Johnson
# setting up the sampler

source('gpriorSetup.R')

w <- c(30, 35, 40, 45, 50)

# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
# for( chain in seq_along(chainSamples) ) {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      # updating g
      temp <- cucumber::slice_sampler_stepping_out(g, log_target,
      w = w[chain], log = TRUE, max = Inf)
      g <- temp$x
      nEval <- temp$nEvaluations
      saving_updates()
    }
    burnin_thinning()
  })
  save_time()
}

# saving samples
saveRDS(chainSamples, file = 'data/steppingOutExploration.rds')

print('finished')
