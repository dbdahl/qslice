args <- commandArgs(trailingOnly = TRUE)
data_use <- args[1]
ii <- as.numeric(args[2])  # job id
dte <- as.numeric(args[3])

##### for testing
# ii <- 54
# dte <- 260518
# data_use <- "mtcars"
# data_use <- "db" # diabetes n = 442
# data_use <- "db40" # diabetes n = 40
# data_use <- "riboflavin"
#####

library("qslice")
library("coda")
source("functions/MH_samplers.R")
source("functions/full_conditionals.R")
source("functions/mcmc_horseshoe.R")
source("functions/tune.R")

load(paste0("schedule_all_", data_use, ".rda"))
source("0_data.R")
source("0_prior.R")

n_iter <- 30e3

run_id <- job_order[ii]
(run_info <- sched[run_id,])
target <- as.character(run_info$target)
(logscale_tau2 <- grepl("log", target))
type <- as.character(run_info$type)
subtype <- as.character(run_info$subtype)

set.seed(dte + ii)

state0 <- list(iter = 0)
state0$sig2 <- prior$s02
state0$beta <- rnorm(dat$p, mean = dat$beta_hat, sd = 0.2)
state0$lam2 <- runif(dat$p, min = 0.0, max = 1.0)
state0$tau2 <- runif(1, min = 0.0, max = 1.0)
state0

sampler0 <- list()

sampler0$tau2 <- list(type = "stepping", subtype = NA, logscale = logscale_tau2,
                     support = c(0.0, 1.0e4),
                     w = ifelse(logscale_tau2, 4.0, 0.005)) # always do initial with a stepping-out slice sampler; dependent on data
## diabetes ltau2 w = 4 always; tau2 w = .003 for full data and .02 (or .003 is just as good...) for n = 40
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

### Additional prep

tau2_samples <- sapply(mc_out$sims, function(x) x$tau2)
# plot(as.mcmc(tau2_samples))
# plot(as.mcmc(log(tau2_samples)))
# coda::effectiveSize(log(tau2_samples))

beta_samples <- sapply(mc_out$sims, function(x) x$beta) |> t()
# indx_check <- 1:4
# (indx_check <- which(colMeans(abs(beta_samples) > 0.01) > 0.5))
(indx_check <- which(apply(beta_samples, 2, function(b) abs(mean(b > 0.0) - 0.5)) > 0.35))
# plot(as.mcmc(beta_samples[,indx_check]))

llam2_samples <- sapply(mc_out$sims, function(x) log(x$lam2)) |> t()
# plot(as.mcmc(llam2_samples[,indx_check]))

# sig2_samples <- sapply(mc_out$sims, function(x) x$sig2)
# plot(as.mcmc(sig2_samples))

# plot(log(tau2_samples), rowMeans(llam2_samples))

# llam2_q95 <- apply(llam2_samples, 1, function(x) quantile(x, 0.95))
# llam2_q95 <- apply(llam2_samples, 1, function(x) mean(x[x > quantile(x, 0.90)]))

# plot(log(tau2_samples), llam2_q95)
# cor(log(tau2_samples), llam2_q95)


if (isTRUE(logscale_tau2)) {
  samples_use <- log(tau2_samples)
} else {
  samples_use <- tau2_samples
}

if (type %in% c("rw", "stepping", "latent")) { # will require tuning

  tune_bnds_init <- c(0.2 * sd(samples_use), 0.8 * diff(range(samples_use)))

  if (type == "latent") {
    state2$latent_s <- mean(tune_bnds_init)
    tune_bnds_init <- rev(1.0 / tune_bnds_init) # uses rate parameter
  }

  sampler_tuned$tau2$type <- type
  sampler_tuned$tau2$subtype <- subtype

  mc_tune <- tune(state = state2, prior = prior, data = dat,
                  sampler = sampler_tuned,
                  param = "tau2",
                  ess_log = TRUE,
                  bnds_init = tune_bnds_init,
                  n_iter = 1e3, n_grid = 3, n_rep = 2,
                  n_rounds = 5, range_frac = 0.67,
                  verbose = TRUE)

  # plot(mc_tune$lvals, mc_tune$lesps)
  # abline(v = log(mc_tune$val_opt), lty = 2)

  sampler_tuned$tau2 <- mc_tune$sampler$tau2

  ### burn-in again with new sampler (so all timing runs begin after iter 10,000)
  mc_out <- mcmc_hs(state = state0, prior = prior, data = dat,
                    sampler = sampler0, # let the slice sampler get to stationarity quickly
                    n_iter = 2e3,
                    save = FALSE, prog = 1000)

  (state1 <- mc_out$state)

  if (type == "latent") {
    state1$latent_s <- mean(tune_bnds_init)
  }

  mc_out <- mcmc_hs(state = state1, prior = prior, data = dat,
                    sampler = sampler_tuned, # let the tuned sampler do the rest
                    n_iter = 8e3,
                    save = FALSE, prog = 1000)

  (state2 <- mc_out$state)

} else if (grepl("samples$", subtype)) { # tune with samples marginally

  (util_type <- substr(subtype, 1, 3))

  tmp_pseu <- pseudo_opt(samples = samples_use,
                         type = "samples",
                         family = "t",
                         degf = c(1, 5),
                         lb = ifelse(logscale_tau2, -Inf, 0.0),
                         ub = Inf,
                         utility_type = util_type,
                         plot = FALSE)

  sampler_tuned$tau2$type <- type
  sampler_tuned$tau2$subtype <- subtype
  sampler_tuned$tau2$pseudo <- tmp_pseu$pseudo
  sampler_tuned$tau2$loc <- tmp_pseu$pseudo$params$loc
  sampler_tuned$tau2$sc <- tmp_pseu$pseudo$params$sc
  sampler_tuned$tau2$degf <- tmp_pseu$pseudo$params$degf
  sampler_tuned$tau2$txt <- tmp_pseu$pseudo$txt

} else if (grepl("samples_reg", subtype)) {

  bqr <- best_quantile_reg(y = log(tau2_samples), X = llam2_samples)

  sampler_tuned$tau2$type <- type
  sampler_tuned$tau2$subtype <- subtype
  sampler_tuned$tau2$bqr <- bqr
  sampler_tuned$tau2$degf <- 1.0

}

### timing run
mc_time <- time_hs(state = state2, prior = prior, data = dat,
                   sampler = sampler_tuned, n_iter = n_iter, param = "tau2",
                   ess_log = TRUE)

mc_time$timing
(samp_p_sec <- mc_time$timing$EffSamp / mc_time$timing$userTime)


### diagnostics
draws_tau2 <- mc_time$draws # are always tau2 (not log scale), but ESS calculation is always for log(tau2)
(tau2_SE <- unname(summary(as.mcmc(draws_tau2))$statistics["Time-series SE"]))
(tau2_mn <- mean(draws_tau2))
(tau2_sd <- sd(draws_tau2))

draws_ltau2 <- log(draws_tau2)
(ltau2_SE <- unname(summary(as.mcmc(draws_ltau2))$statistics["Time-series SE"]))
(ltau2_mn <- mean(draws_ltau2))
(ltau2_sd <- sd(draws_ltau2))

# plot(as.mcmc(draws_tau2))
# plot(as.mcmc(draws_ltau2))

if (type %in% c("rw", "stepping", "latent")) {
  tune_param <- mc_tune$val_opt
} else {
  if (subtype == "samples_reg") {
    tune_param <- bqr$qntle
  } else {
    tune_param <- NA
  }
}

if (type == "Qslice") {
  n_iter_qs <- 2e3
  mc_out <- mcmc_hs(state = state2, prior = prior, data = dat,
                        sampler = sampler_tuned,
                        n_iter = n_iter_qs,
                        save = TRUE, prog = 0)

  draws_u <- sapply(mc_out$extras, function(x) x$u)
  (IAT_u <- unname(n_iter_qs / effectiveSize(draws_u)))
  (AUC <- auc(u = draws_u))
} else {
  IAT_u <- NA
  AUC <- NA
}

tempDf <- data.frame(target = target,
                     type = type,
                     subtype = subtype,
                     rep = run_info$rep,
                     run_id = run_id,
                     ii = ii,
                     n_iter = n_iter,
                     tau2_mn = tau2_mn,
                     tau2_sd = tau2_sd,
                     tau2_SE = tau2_SE,
                     ltau2_mn = ltau2_mn,
                     ltau2_sd = ltau2_sd,
                     ltau2_SE = ltau2_SE,
                     nEval = mc_time$timing$nEval,
                     ESS = mc_time$timing$EffSamp,
                     userTime = mc_time$timing$userTime,
                     sysTime = mc_time$timing$sysTime,
                     elapsedTime = mc_time$timing$elapsedTime,
                     sampPsec = samp_p_sec,
                     tuneParam = tune_param,
                     iat_u = IAT_u,
                     auc = AUC)

write.table(tempDf, file = paste0("output/",
                                  "target", target,
                                  "_type", type,
                                  "_subtype", subtype,
                                  "_rep", run_info$rep,
                                  "_dte", dte,
                                  ".csv"),
            append = FALSE, sep = ",",
            row.names = FALSE, col.names = TRUE)

print('finished')


quit(save = "no")
