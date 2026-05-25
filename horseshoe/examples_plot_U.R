rm(list=ls()); dev.off()

set.seed(1)

library("qslice")
library("coda")
source("functions/MH_samplers.R")
source("functions/full_conditionals.R")
source("functions/mcmc_horseshoe.R")
source("functions/tune.R")

# data_use <- "db40"
data_use <- "db"
source("0_data.R")
source("0_prior.R")

state0 <- list(iter = 0)
state0$sig2 <- prior$s02
state0$beta <- rnorm(dat$p, mean = dat$beta_hat, sd = 0.2)
state0$lam2 <- runif(dat$p, min = 0.0, max = 1.0)
state0$tau2 <- runif(1, min = 0.0, max = 1.0)
state0

n_iter <- 20e3


logscale_tau2 <- TRUE


sampler <- list()


### AUC Samples

## gather preliminary samples
sampler0 <- list()

sampler0$tau2 <- list(type = "stepping", subtype = NA, logscale = logscale_tau2,
                      support = c(0.0, 1.0e4),
                      w = ifelse(logscale_tau2, 4.0, 0.005)) # always do initial with a stepping-out slice sampler; dependent on data
sampler_tuned <- sampler0 # make a copy

### initial burn-in
mc_out <- mcmc_hs(state = state0, prior = prior, data = dat,
                  sampler = sampler0,
                  n_iter = 2e3,
                  save = FALSE, prog = 500)

(state1 <- mc_out$state)

### continue the chain and save samples for tuning
mc_out <- mcmc_hs(state = state1, prior = prior, data = dat,
                  sampler = sampler0,
                  n_iter = 2e3, n_thin = 4,
                  save = TRUE, prog = 500) # samples for pseudo-target

(state2 <- mc_out$state)

tau2_samples <- sapply(mc_out$sims, function(x) x$tau2)
llam2_samples <- sapply(mc_out$sims, function(x) log(x$lam2)) |> t()

llam2_q95 <- apply(llam2_samples, 1, function(x) quantile(x, 0.95))

pdf(file = paste0("plots/ltau2_llam2_reg_", data_use, ".pdf"), width = 1.2*4, height = 4)
par(mar = c(4.3, 4.3, 1.1, 1.1))
plot(llam2_q95, log(tau2_samples),
     ylab = expression(log(tau^2)),
     xlab = bquote("95th percentile of"~ log(lambda^2)),
     axes = FALSE, cex = 0.15)
axis(side = 1)
axis(side = 2)
dev.off()

cor(llam2_q95, log(tau2_samples))

if (isTRUE(logscale_tau2)) {
  samples_use <- log(tau2_samples)
} else {
  samples_use <- tau2_samples
}


## AUC samples

util_type <- "AUC"

tmp_pseu <- pseudo_opt(samples = samples_use,
                       type = "samples",
                       family = "t",
                       degf = c(1, 5),
                       lb = ifelse(logscale_tau2, -Inf, 0.0),
                       ub = Inf,
                       utility_type = util_type,
                       plot = FALSE)

sampler_tuned$tau2$type <- "Qslice"
sampler_tuned$tau2$subtype <- "AUC_samples"
sampler_tuned$tau2$pseudo <- tmp_pseu$pseudo
sampler_tuned$tau2$loc <- tmp_pseu$pseudo$params$loc
sampler_tuned$tau2$sc <- tmp_pseu$pseudo$params$sc
sampler_tuned$tau2$degf <- tmp_pseu$pseudo$params$degf
(sampler_tuned$tau2$txt <- tmp_pseu$pseudo$txt)

mc_out_auc <- mcmc_hs(state = state2, prior = prior, data = dat,
                  sampler = sampler_tuned,
                  n_iter = n_iter,
                  save = TRUE, prog = 1000)

draws_u_auc <- sapply(mc_out_auc$extras, function(x) x$u)

gthm <- theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11))

ggplot(data.frame(u = draws_u_auc), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_classic() + gthm +
  ggtitle(label = paste0("Method: AUC - samples\n", "AUC: ", round(auc(u = draws_u_auc), 2),
                         "     IAT: ", round(n_iter / effectiveSize(draws_u_auc), 2)))

ggsave(file = paste0("plots/histU_AUCsamples_", data_use, ".pdf"), height = 2.5, width = 3)


## MSW samples

util_type <- "MSW"

tmp_pseu <- pseudo_opt(samples = samples_use,
                       type = "samples",
                       family = "t",
                       degf = c(1, 5),
                       lb = ifelse(logscale_tau2, -Inf, 0.0),
                       ub = Inf,
                       utility_type = util_type,
                       plot = FALSE)

sampler_tuned$tau2$type <- "Qslice"
sampler_tuned$tau2$subtype <- "MSW_samples"
sampler_tuned$tau2$pseudo <- tmp_pseu$pseudo
sampler_tuned$tau2$loc <- tmp_pseu$pseudo$params$loc
sampler_tuned$tau2$sc <- tmp_pseu$pseudo$params$sc
sampler_tuned$tau2$degf <- tmp_pseu$pseudo$params$degf
(sampler_tuned$tau2$txt <- tmp_pseu$pseudo$txt)

mc_out_msw <- mcmc_hs(state = state2, prior = prior, data = dat,
                  sampler = sampler_tuned,
                  n_iter = n_iter,
                  save = TRUE, prog = 1000)

draws_u_msw <- sapply(mc_out_msw$extras, function(x) x$u)

gthm <- theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11))

ggplot(data.frame(u = draws_u_msw), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_classic() + gthm +
  ggtitle(label = paste0("Method: MSW - samples\n", "AUC: ", round(auc(u = draws_u_msw), 2),
                         "     IAT: ", round(n_iter / effectiveSize(draws_u_msw), 2)))

ggsave(file = paste0("plots/histU_MSWsamples_", data_use, ".pdf"), height = 2.5, width = 3)


## Reg samples

bqr <- best_quantile_reg(y = log(tau2_samples), X = llam2_samples)

sampler_tuned$tau2$type <- "Qslice"
sampler_tuned$tau2$subtype <- "samples_reg"
sampler_tuned$tau2$bqr <- bqr
sampler_tuned$tau2$pseu_bqr <- pseudo_opt(samples = residuals(bqr$model),
                                          type = "samples", famil = "t",
                                          degf = c(1, 5), utility_type = "AUC",
                                          plot = FALSE)

mc_out_reg <- mcmc_hs(state = state2, prior = prior, data = dat,
                  sampler = sampler_tuned,
                  n_iter = n_iter,
                  save = TRUE, prog = 1000)

draws_u_reg <- sapply(mc_out_reg$extras, function(x) x$u)

gthm <- theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11))

ggplot(data.frame(u = draws_u_reg), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  xlab(expression(psi)) + theme_classic() + gthm +
  ggtitle(label = paste0("Method: Regression\n", "AUC: ", round(auc(u = draws_u_reg), 2),
                         "     IAT: ", round(n_iter / effectiveSize(draws_u_reg), 2)))

ggsave(file = paste0("plots/histU_REGsamples_", data_use, ".pdf"), height = 2.5, width = 3)
