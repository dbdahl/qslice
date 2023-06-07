## g prior sampling using transform slice psuedo updated at every iteration using optim based on samples
## Author: Sam Johnson

source('gpriorSetup.R')

# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      ##
      if( i == 1 ) {
        approxSamples <- approx_samples(x = g, lf = log_target, samples = 5000, w = 25)
        psuedoFit <- opt_Cauchy_auc_data(approxSamples, lb = 0)
        psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
      }
      if(i %% 1 == 1000 & i > 2) {
        approxSamples <- chainSamples[[chain]]$g[1:(i-1)]
        psuedoFit <- opt_Cauchy_auc_data(approxSamples, lb = 0)
        psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
      }
      ##
      temp <- cucumber::slice_sampler_transform(x = g, target = log_target,
                                                pseudo_log_pdf = psuedoTarget$pseudo_log_pdf,
                                                pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)
      g <- temp$
      nEval <- temp$nEvaluations
      u <- temp$u
      saving_updates(method = 'Transform')
    }
    burnin_thinning(method = 'Transform')
  })
  save_time()
}

# saving samples
saveRDS(chainSamples, file = 'data/transformEveryIterOptimSamples.rds')


print('finished')

