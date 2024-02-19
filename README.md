# Cucumber Package

Functions to facilitate implementing the quantile slice sampler and other popular slice samplers. Package includes utility functions for specifying psueudo-targets, for diagnostics, and for tuning.

This repository has two parts found in two different directories. The **cucumber** directory has the code needed to run each slice sampler (GESS, Latent, Stepping Out, Quantile). The **samCode** directory has all the code to run the simulations used to compare each of the sampling methods.

The **samCode** directory contains sub-directories for the G prior simulation study and each of the standard targets that were tested, standard normal, gamma, and inverse gamma. Each sub-directory has README files that explain each script.

The following code gives a demonstration of the typical workflow. Note that `r slice_sampler_transform()` can be embedded within any Gibbs sampler.

``` r
## use gamma(shape = alpha) as a target
alpha <- 2.5
target <- list(d = function(x) x^(alpha - 1.0) * exp(-x) * (x > 0), # unnormalized density
               ld = function(x) (alpha - 1.0)*log(x) - x + ifelse(x > 0.0, 0.0, -Inf) # log unnormalized density
               )

## find pseudo-target that maximizes AUC; plots the transformed target
pseudo <- opt_t(target, type = "function", lb = 0.0, 
                use_meanSliceWidth = FALSE,
                verbose = FALSE, plot = TRUE)
pseudo$pseu$t
pseudo$util

## compare pseudo-target and target densities (note target is not normalized)
curve(target$d(x), from = -1, to = 4*alpha, n = 1e3)
curve(pseudo$pseu$d(x), from = -1, to = 4*alpha, add = TRUE, lty = 2, n = 1e3)
legend("topright", lty = 1:2, legend = c("target", "pseudo-target"), bty = "n")

## set up MCMC
n_sim <- 10e3
samp_x <- samp_psi <- numeric(n_sim)
samp_x[1] <- 0.5  # initialize
samp_psi[1] <- pseudo$pseu$p(samp_x[1])
n_eval <- 0  # count target evaluations

## run quantile slice sampler
for (i in 2:n_sim) {
  state <- slice_sampler_transform(samp_x[i-1], 
                                   target = target$ld, 
                                   pseudo_log_pdf = pseudo$pseu$ld, 
                                   pseudo_inv_cdf = pseudo$pseu$q,
                                   log = TRUE)
  n_eval <- n_eval + state$nEvaluations
  samp_psi[i] <- state$u
  samp_x[i] <- state$x
}

## check samples
n_eval / (n_sim - 1)  # target evaluations per iteration of MCMC
hist(samp_x, freq = FALSE, n = 20)
curve(dgamma(x, alpha), col = "blue", lwd = 2, add = TRUE)
ks.test(samp_x, pgamma, shape = alpha)  # null hypothesis: samp_x ~iid gamma(alpha)

## diagnostics
hist(samp_psi)  # want close to uniform
utility_shrinkslice(u = samp_psi, type = "samples")  # want AUC close to 1

## create a new pseudo-target based on first round of samples (tune)
pseudo2 <- opt_t(samples = samp_x, type = "samples", lb = 0.0, 
                use_meanSliceWidth = FALSE,
                verbose = FALSE, plot = TRUE)
pseudo2$pseu$t
```
