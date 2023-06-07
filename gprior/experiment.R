## expermiment


source('gpriorSetup.R')


# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(Nsamples)) {
      update_parameters()
      ## RAND WALK
      # temp <- randWalk(int.x = g, lf = log_target, c = 30)
      # nEval <- temp$accept
      ## STEPPING OUT
      # temp <- cucumber::slice_sampler_stepping_out(g, log_target,
      #                                              w = 50, log = TRUE, max = Inf)
      ## LATENT
      # temp <- cucumber::slice_sampler_latent(x = g, s = 1, log_target,
      #                                                rate = 0.001, log = TRUE)
      ## GESS
      # mu = 20, sigma = 15, df = 3
      # temp <- cucumber::slice_sampler_generalized_elliptical(x = g, log_target,
      #                                                        mu = 20, sigma = 300, df = 3, log = TRUE)
      ## TRANSFORM
      # LAPLACE
      # if(i == 1) {
        # psuedoFit <- lapproxt(target, 10, lb = 0, maxub = 200)
        # psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = 0, name = 'Laprox')
      # }
      # temp <- cucumber::slice_sampler_transform(x = g, target = log_target,
      #                                           pseudo_log_pdf = psuedoTarget$pseudo_log_pdf,
      #                                           pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)
      # AUTO
      # if(i == 1) {
      #   approxSamples <- approx_samples(x = g, lf = log_target, samples = 5000, w = 25)
      #   psuedoFit <- fit_trunc_Cauchy(approxSamples, lb = 0)
      #   psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
      # }
      # temp <- cucumber::slice_sampler_transform(x = g, target = log_target,
      #                                           pseudo_log_pdf = psuedoTarget$pseudo_log_pdf,
      #                                           pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)
      # # OPTIM SAMPLES
      if(i == 1) {
        approxSamples <- approx_samples(x = g, lf = log_target, samples = 5000, w = 25)
        psuedoFit <- opt_Cauchy_auc_data(approxSamples, lb = 0)
        psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
      }
      temp <- cucumber::slice_sampler_transform(x = g, target = log_target,
                                                pseudo_log_pdf = psuedoTarget$pseudo_log_pdf,
                                                pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)
      g <- temp$x
      nEval <- temp$nEvaluations
      
      saving_updates()
    }
    burnin_thinning()
  })
  save_time()
}

gSamples <- sapply(chainSamples, FUN = \(list) list$g)

traceplot(gSamples)

resultsFunc(chainSamples = chainSamples)

# 

# Laplace
psuedoFit <- lapproxt(target, 10, lb = 0, maxub = 200)
psuedoTarget1 <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = 0, name = 'Laprox')
# Auto
approxSamples <- approx_samples(x = g, lf = log_target, samples = 5000, w = 25)
psuedoFit <- fit_trunc_Cauchy(approxSamples, lb = 0)
psuedoTarget2 <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
# Optim Samples
psuedoFit <- opt_Cauchy_auc_data(approxSamples, lb = 0)
psuedoTarget3 <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
# xx <- seq(from = 0, to = 200, length.out = 1000)
# yy <- target(xx)
# yy1 <- sapply(xx, FUN = \(x) exp(psuedoTarget1$pseudo_log_pdf(x)))
# yy2 <- sapply(xx, FUN = \(x) exp(psuedoTarget2$pseudo_log_pdf(x)))
# yy3 <- sapply(xx, FUN = \(x) exp(psuedoTarget3$pseudo_log_pdf(x)))
# 
# plot(xx,yy, xlim = c(0,600), ylim = c(0,9e-3))
# points(xx,yy1, col = 'red')
# points(xx,yy2, col = 'darkgreen')
# points(xx, yy3, col = 'dodgerblue4')
# 
# 
# ## group 6
# xx = seq(0, 200, length=1001)
# uu = seq(0, 1, length=1001)
# truth = list(d = function(x) {target(x)},
#              q = function(u) {qt(u, df=5)},
#              t = "target")
# 
# pseu = list()
# pseu[[1]] = list(d = function(x) {exp(psuedoTarget1$pseudo_log_pdf(x))},
#                  q = function(u) {psuedoTarget1$pseudo_inv_cdf(u)},
#                  t = "Laplace")
# 
# pseu[[2]] = list(d = function(x) {exp(psuedoTarget2$pseudo_log_pdf(x))},
#                  q = function(u) {psuedoTarget2$pseudo_inv_cdf(u)},
#                  t = "Auto")
# 
# pseu[[3]] = list(d = function(x) {exp(psuedoTarget3$pseudo_log_pdf(x))},
#                  q = function(u) {psuedoTarget3$pseudo_inv_cdf(u)},
#                  t = "Optim Samples")
# 
# ### plot
# (n_pseu = length(pseu))
# 
# m = matrix(c(1,2,3), nrow=1)
# layout(mat=m, widths=c(1, 1, 0.7))
# par(mar=c(4,1,1,1))
# 
# plot(xx, truth$d(xx), type="l", lwd=3,
#      ylim=c(0,9e-3),
#      xlab=expression(theta), axes = FALSE, col = 'black')
# axis(side=1)
# for(i in 1:n_pseu) {
#   lines(xx, pseu[[i]]$d(xx), col=i+1, lwd=3, lty=i+1)	
# }
# 
# plot(uu, rep(1.0, length(uu)), type="l", col=1, lwd=3, 
#      ylim=c(0,0.1),
#      xlab=expression(psi), ylab="", axes=FALSE, bty="L")
# axis(side=1)
# 
# for(i in 1:n_pseu) {
#   lines(uu, truth$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col=i+1, lwd=2, lty=i+1)
# }
# 
# plot(1, type="n", axes=FALSE, xlab="", ylab="")
# legend("left", bty="n", inset=c(0.0, 0.0),
#        legend=c(truth$t, sapply(pseu, function(x) x$t)),
#        col=1:(n_pseu+1), lwd=2, lty=1:(n_pseu+1), cex = 1.2)
# psuedoTarget1
# psuedoTarget2
# psuedoTarget3
# 
