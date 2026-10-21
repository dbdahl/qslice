# Qslice Package

Functions to facilitate implementing the quantile slice sampler and other popular slice samplers. Package includes utility functions for specifying pseudo-targets, for diagnostics, and for tuning.

The **qslice** directory contains the R package, which implements the quantile slice sampler and related samplers samplers.

The other directories hold earlier development code can largely be ignored.  Code to reproduce the simulation studies in the paper is in the [qslice_examples](https://github.com/mheiner/qslice_examples) repository.

The following code gives a demonstration of the typical workflow. Note that `slice_quantile()` can be embedded within any Gibbs sampler.

``` r
## use gamma(shape = alpha) as a target
alpha <- 2.5
target <- list(d = function(x) ifelse(x > 0.0, x^(alpha - 1.0) * exp(-x), 0.0), # unnormalized density
               ld = function(x) ifelse(x > 0.0, (alpha - 1.0)*log(x) - x, -Inf) # log unnormalized density
               )

## find pseudo-target that maximizes AUC; plots the transformed target
par(mfrow = c(1, 2))
pseudo <- pseudo_opt(log_target = target$ld, 
                     type = "function", 
                     family = "t",
                     lb = 0.0, 
                     utility_type = "AUC",
                     verbose = FALSE, plot = TRUE)
pseudo$pseudo$txt
pseudo$utility

## set up MCMC
n_iter <- 10e3
samp_x <- samp_psi <- numeric(n_iter + 1)
samp_x[1] <- 0.5  # initialize
samp_psi[1] <- pseudo$pseudo$p(samp_x[1])
n_eval <- 0  # count target evaluations

## run quantile slice sampler
for (i in 2:(n_iter+1)) {
  state <- slice_quantile(samp_x[i-1], 
                          log_target = target$ld, 
                          pseudo = pseudo$pseudo)
  n_eval <- n_eval + state$nEvaluations
  samp_psi[i] <- state$u
  samp_x[i] <- state$x
}

## check samples
n_eval / n_iter  # target evaluations per iteration of MCMC
hist(samp_x, freq = FALSE, n = 20)
curve(dgamma(x, alpha), col = "blue", lwd = 2, add = TRUE)
ks.test(samp_x, pgamma, shape = alpha)  # null hypothesis: samp_x ~iid gamma(alpha)

## diagnostics
hist(samp_psi)  # want close to uniform
utility_shrinkslice(u = samp_psi, utility_type = "AUC")  # want AUC close to 1

## create a new pseudo-target based on first round of samples (tune)
pseudo2 <- pseudo_opt(samples = samp_x, 
                      type = "samples",
                      family = "t",
                      lb = 0.0, 
                      utility_type = "AUC",
                      verbose = FALSE, plot = TRUE)
pseudo2$pseudo$txt
```
