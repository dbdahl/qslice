## g prior sampling using transform slice just one psuedo prior fuond using optim
## Author: Sam Johnson

source('gpriorSetup.R')

# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      if(i == 1) {
        truth = list(ld = target,
                     t = "Full Conditional of g",
                     lb = 0,
                     ub = 200)
        truth$dld = function(x) numDeriv::grad(truth$ld, x=x)
        psuedoFit <- opt_Cauchy_auc(truth)
        psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$par[1], sc = psuedoFit$par[2], lb = 0, ub = Inf)
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
saveRDS(chainSamples, file = 'data/transformOnceOptim.rds')

print('finished')

