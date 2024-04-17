rm(list=ls()); dev.off()

set.seed(1)

library("cucumber")
source("0_data.R")
source("0_prior.R")
source("functions/MH_samplers.R")
source("functions/full_conditionals.R")
source("functions/mcmc_gprior.R")
source("functions/tune.R")
source("functions/Laplace_g.R")
source("functions/invgamma.R")

state <- list(iter = 0)
state$psi <- prior$a_0 / prior$b_0
state$g <- prior$g_max / 2
state$beta <- prior$beta_0
state

n_iter <- 10e3

sampler <- list()

sampler$g <- list(type = "Gibbs", subtype = NA, logG = F)

sampler$g <- list(type = "rw", subtype = NA, logG = F, c = 65.0) # 30-3600; 20-2200; 25-3000; 35-4500; 40-4800; 45-5300; 50-5500; 55-5600; 60-5700; 65-5750; 70-5600; 75-5400
sampler$g <- list(type = "stepping", subtype = NA, logG = F, w = 70.0) # 30-13000; 20-12800; 40-13500; 45-13600; 50-13700; 60-13600; 70-13500 volatile
sampler$g <- list(type = "gess", subtype = NA, logG = F, loc = 40.0, sc = 30, degf = 3)
sampler$g <- list(type = "latent", subtype = NA, logG = F, rate = 0.001); state$latent_s <- 5.0 # bnds_init = c(0.001, 0.2); 4 rounds ok

sampler$pseu <- pseudo_t_list(loc = 30, sc = 20, degf = 1, lb = 0.0, ub = prior$g_max)
sampler$g <- list(type = "imh", subtype = NA, logG = F,
                 pseudo_lpdf = sampler$g$pseudo_lpdf,
                 # pseudo_lpdf = sampler$pseu$ld,
                 pseudo_inv_cdf = sampler$g$pseudo_inv_cdf)
                 # pseudo_inv_cdf = sampler$pseu$q)
sampler$g <- list(type = "Qslice", subtype = "manual",
                  logG = F,
                 pseudo_lpdf = sampler$pseu$ld,
                 pseudo_inv_cdf = sampler$pseu$q)

## I think we'll drop these
sampler$g <- list(type = "Qslice", subtype = "Laplace", sc_adj = 1.0, logG = F)
sampler$g <- list(type = "Qslice", subtype = "Laplace_wide", sc_adj = 1.5, maxit = 2500, logG = F)
##

sampler$g <- list(type = "Qslice", subtype = "Laplace_analytic", sc_adj = 1.0, degf = 1, logG = F)
sampler$g <- list(type = "Qslice", subtype = "Laplace_analytic_wide", sc_adj = 1.5, degf = 1, logG = F)

sampler$g <- list(type = "imh", subtype = "Laplace_analytic", sc_adj = 1.0, degf = 1, logG = F)
sampler$g <- list(type = "imh", subtype = "Laplace_analytic_wide", sc_adj = 1.5, degf = 1, logG = F)

sampler$g <- list(type = "Qslice", subtype = "MM", sc_adj = 1.0, degf = 1, logG = F)
sampler$g <- list(type = "Qslice", subtype = "MM_wide", sc_adj = 1.5, degf = 1, logG = F)

sampler$g <- list(type = "imh", subtype = "MM", sc_adj = 1.0, degf = 1, logG = F)
sampler$g <- list(type = "imh", subtype = "MM_wide", sc_adj = 1.5, degf = 1, logG = F)



sampler$g <- list(type = "rw", subtype = NA, logG = T, c = 1.5) # bnds_init = c(0.5, 5.0)
sampler$g <- list(type = "stepping", subtype = NA, logG = T, w = 3.0) # bnds_init = c(0.5, 5.0)
sampler$g <- list(type = "gess", subtype = NA, logG = T, loc = 3.4, sc = 0.6, degf = 5)
sampler$g <- list(type = "latent", subtype = NA, logG = T, rate = 0.01); state$latent_s <- 5.0  # bnds_init = c(0.001, 0.2); 4 rounds ok

sampler$pseu <- pseudo_t_list(loc = 3.4, sc = 0.6, degf = 5, lb = 0.0, ub = prior$g_max)
sampler$g <- list(type = "imh", subtype = NA, logG = T,
                 pseudo_lpdf = sampler$pseu$ld,
                 pseudo_inv_cdf = sampler$pseu$q)
sampler$g <- list(type = "Qslice", logG = T,
                 pseudo_lpdf = sampler$pseu$ld,
                 pseudo_inv_cdf = sampler$pseu$q)


### Sampling
mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = 2e3,
                      save = TRUE, prog = 0)

(state <- mc_out$state)

mc_tune <- tune(state = state, prior = prior, data = dat,
                sampler = sampler,
                bnds_init = c(5.0, 100.0),
                n_iter = 1000, n_grid = 5, n_rep = 3,
                n_rounds = 5, range_frac = 0.5,
                verbose = TRUE)

(sampler$g <- mc_tune$sampler$g)

tmp_pseu <- opt_t(samples = sapply(mc_out$sims, function(x) x$g),
                  type = "samples",
                  lb = 0.0,
                  ub = prior$g_max,
                  coeffs = c(1, 0),
                  degf = c(1, 5),
                  plot = TRUE)

sampler$g <- list(type = "Qslice",
                  subtype = "AUC_samples",
                  logG = sampler$g$logG,
                  pseudo_lpdf = tmp_pseu$pseu$ld,
                  pseudo_inv_cdf = tmp_pseu$pseu$q,
                  t = tmp_pseu$pseu$t)

mc_time <- time_gprior(state = state, prior = prior, data = dat,
                       sampler = sampler, n_iter = n_iter)
# hist(mc_time$draws, breaks = 50)
mc_time$timing
mc_time$timing$EffSamp / mc_time$timing$userTime

library("coda")
library("ggplot2")

names(mc_out)
draws_g <- sapply(mc_out$sims, function(x) x$g)
plot(as.mcmc(draws_g))
hist(draws_g, breaks = 50, freq = FALSE)

## diagnostics
draws_u <- sapply(mc_out$extras, function(x) x$u)
hist(draws_u, breaks = 30)

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_bw() + ggtitle(label = paste0("Laplace analytic - wide\n", "AUC: ", round(auc(u = draws_u), 2)))

ggsave(file = "plots/histU_LaplaceAnalyticWide.pdf", height = 2.5, width = 3)

hist(draws_g, breaks = 50, freq = FALSE); curve(exp(sampler$g$pseudo_lpdf(x)),
                                                from = 0, to = prior$g_max,
                                                col = "red", add = TRUE)


sapply(mc_out$extras, function(x) x$accept) |> mean()


draws_psi <- sapply(mc_out$sims, function(x) x$psi)
plot(as.mcmc(draws_psi))

draws_beta <- sapply(mc_out$sims, function(x) x$beta) |> t()
colnames(draws_beta) <- colnames(dat$X)
plot(as.mcmc(draws_beta))
summary(draws_beta)
colMeans(draws_beta)
dat$beta_mle
