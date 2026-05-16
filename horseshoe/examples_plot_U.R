rm(list=ls()); dev.off()

set.seed(1)

library("qslice")
library("coda")
library("ggplot2")

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


### AUC Samples

## gather preliminary samples
sampler$g <- list(type = "stepping", subtype = NA, logG = F, w = 70.0)

mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = 2e3,
                      save = TRUE, prog = 0)

(state <- mc_out$state)

tmp_pseu <- pseudo_opt(samples = sapply(mc_out$sims, function(x) x$g),
                       type = "samples",
                       family = "t",
                       degf = c(1, 5),
                       lb = 0.0,
                       ub = prior$g_max,
                       plot = TRUE,
                       nbins = 20)

sampler$g <- list(type = "Qslice",
                  subtype = "AUC_samples",
                  logG = sampler$g$logG,
                  pseudo = tmp_pseu$pseudo,
                  t = tmp_pseu$pseu$t)

mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = n_iter,
                      save = TRUE, prog = 0)

(state <- mc_out$state)

draws_u <- sapply(mc_out$extras, function(x) x$u)

gthm <- theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11))

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_classic() + gthm +
  ggtitle(label = paste0("AUC - samples\n", "AUC: ", round(auc(u = draws_u), 2)))

ggsave(file = "plots/histU_AUCsamples.pdf", height = 2.5, width = 3)



### Laplace Analytic

sampler$g <- list(type = "Qslice", subtype = "Laplace_analytic", sc_adj = 1.0, degf = 1, logG = F)

mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = n_iter,
                      save = TRUE, prog = 0)

(state <- mc_out$state)

draws_u <- sapply(mc_out$extras, function(x) x$u)

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_classic() + gthm +
  ggtitle(label = paste0("Laplace analytic\n", "AUC: ", round(auc(u = draws_u), 2)))

ggsave(file = "plots/histU_LaplaceAnalytic.pdf", height = 2.5, width = 3)



### Laplace Analytic Wide

sampler$g <- list(type = "Qslice", subtype = "Laplace_analytic_wide", sc_adj = 1.5, degf = 1, logG = F)

mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = n_iter,
                      save = TRUE, prog = 0)

(state <- mc_out$state)

draws_u <- sapply(mc_out$extras, function(x) x$u)

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_classic() + gthm +
  ggtitle(label = paste0("Laplace analytic - wide\n", "AUC: ", round(auc(u = draws_u), 2)))

ggsave(file = "plots/histU_LaplaceAnalyticWide.pdf", height = 2.5, width = 3)

