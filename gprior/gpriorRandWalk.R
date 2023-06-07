## g prior sampling using rand walk
## Author: Sam Johnson

# setting up the sampler
source('gpriorSetup.R')


# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      temp <- randWalk(int.x = g, lf = log_target, c = 30)
      g <- temp$x
      nEval <- temp$accept
      saving_updates()
    }
    burnin_thinning()
  })   
  save_time()
}

# finding the acceptance rate
acceptanceRate <- c(sapply(chainSamples, FUN = \(list) list$nEval)) |> mean() |> round(2)
print(paste0('The acceptance rate is:', acceptanceRate))

# saving samples
saveRDS(chainSamples, file = 'data/randWalk.rds')

print('finished')
