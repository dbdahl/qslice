## g prior sampling using transform slice psuedo updated once using laplace
## Author: Sam Johnson

source('gpriorSetup.R')

# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
     update_parameters()
      if(i == 1) {
        psuedoFit <- lapproxt(target, 10, lb = 0, maxub = 200)
        psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = 0, name = 'Laprox')
      }
      temp <- cucumber::slice_sampler_transform(x = g, target = log_target,
                                                pseudo_log_pdf = psuedoTarget$pseudo_log_pdf,
                                                pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)
      g <- temp$x
      nEval <- temp$nEvaluations
      u <- temp$u
      saving_updates(method = 'Transform')
    }
    burnin_thinning(method = 'Transform')
  })
  save_time()
}

# saving samples
saveRDS(chainSamples, file = 'data/transformOnceLaplace.rds')


print('finished')
