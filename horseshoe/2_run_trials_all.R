args <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(args[1])  # job id
dte <- as.numeric(args[2])

##### for testing
ii <- 2
dte <- 260511
# data_use <- "mtcars"
data_use <- "diabetes"
# data_use <- "riboflavin"
#####

library("qslice")
library("coda")
source("0_data.R")
source("0_prior.R")
source("functions/MH_samplers.R")
source("functions/full_conditionals.R")
source("functions/mcmc_horseshoe.R")
source("functions/tune.R")

load("schedule_all.rda")

n_iter <- 10e3

run_id <- job_order[ii]
(run_info <- sched[run_id,])
target <- as.character(run_info$target)
(logscale_tau2 <- grepl("log", target))
type <- as.character(run_info$type)
subtype <- as.character(run_info$subtype)

set.seed(dte + ii)

state <- list(iter = 0)
state$sig2 <- prior$s02
state$beta <- rnorm(dat$p, mean = dat$beta_hat, sd = 0.2)
state$lam2 <- runif(dat$p, min = 0.0, max = 1.0)
state$tau2 <- runif(1, min = 0.0, max = 1.0)
state

sampler <- list()

sampler$tau2 <- list(type = "stepping", subtype = NA, logscale = logscale_tau2,
                     w = ifelse(logscale_tau2, 1.0, 0.5)) # always do initial with a stepping-out slice sampler; CHANGE THIS BASED ON data set


### initial burn-in
mc_out <- mcmc_hs(state = state, prior = prior, data = dat,
                  sampler = sampler,
                  n_iter = 1e3,
                  save = FALSE, prog = 100,
                  upper_tau2 = 1.0e9)
                  # appx_threshold = 0.0)

(state <- mc_out$state)

### continue the chain and save samples
mc_out <- mcmc_hs(state = state, prior = prior, data = dat,
                  sampler = sampler,
                  # n_iter = 5e3, n_thin = 5,
                  n_iter = 1e3, n_thin = 5,
                  save = TRUE, prog = 100,
                  upper_tau2 = 1.0e9) # samples for pseudo-target
                  # appx_threshold = 0.0) # samples for pseudo-target

(state <- mc_out$state)

### Additional prep

tau2_samples <- sapply(mc_out$sims, function(x) x$tau2)
plot(as.mcmc(tau2_samples))
plot(as.mcmc(log(tau2_samples)))
coda::effectiveSize(log(tau2_samples))

beta_samples <- sapply(mc_out$sims, function(x) x$beta) |> t()
indx_check <- 1:4
(indx_check <- which(colMeans(abs(beta_samples) > 0.05) > 0.25))
plot(as.mcmc(beta_samples[,indx_check]))

llam2_samples <- sapply(mc_out$sims, function(x) log(x$lam2)) |> t()
plot(as.mcmc(llam2_samples[,indx_check]))

sig2_samples <- sapply(mc_out$sims, function(x) x$sig2)
plot(as.mcmc(sig2_samples))


if (isTRUE(logscale_tau2)) {
  samples_use <- log(tau2_samples)
} else {
  samples_use <- tau2_samples
}

if (type %in% c("rw", "stepping", "latent")) { # will require tuning

  tune_bnds_init <- c(0.2 * sd(samples_use), 0.8*diff(range(samples_use)))

  if (type == "latent") {
    state$latent_s <- mean(tune_bnds_init)
    tune_bnds_init <- rev(1.0 / tune_bnds_init) # uses rate parameter
  }

  sampler$tau2$type <- type
  sampler$tau2$subtype <- subtype

  mc_tune <- tune(state = state, prior = prior, data = dat,
                  sampler = sampler,
                  param = "tau2",
                  bnds_init = tune_bnds_init,
                  n_iter = 1e3, n_grid = 5, n_rep = 3,
                  n_rounds = 5, range_frac = 0.5,
                  verbose = TRUE)

  sampler$tau2 <- mc_tune$sampler$tau2

  if (type == "latent") {
    state <- mc_tune$state # since latent_s was burned in during tuning
  }

} else if (grepl("samples", subtype)) { # tune with samples

  (util_type <- substr(subtype, 1, 3))

  tmp_pseu <- pseudo_opt(samples = samples_use,
                         type = "samples",
                         family = "t",
                         degf = c(1, 5),
                         lb = ifelse(logscale_tau2, -Inf, 0.0),
                         ub = Inf,
                         utility_type = util_type,
                         plot = FALSE)

  sampler$tau2$type <- type
  sampler$tau2$subtype <- subtype
  sampler$tau2$pseudo <- tmp_pseu$pseudo
  sampler$tau2$loc <- tmp_pseu$pseudo$params$loc
  sampler$tau2$sc <- tmp_pseu$pseudo$params$sc
  sampler$tau2$degf <- tmp_pseu$pseudo$params$degf
  sampler$tau2$txt <- tmp_pseu$pseudo$txt

} else if (grepl("Laplace_analytic", subtype)) {

  sampler$tau2$type <- type
  sampler$tau2$subtype <- subtype

  if (isTRUE(logscale_tau2)) {

    sampler$tau2$degf = 5

    if (grepl("wide", subtype)) {
      sampler$tau2$sc_adj <- 1.2
    } else {
      sampler$tau2$sc_adj <- 1.0
    }

  } else {

    sampler$tau2$degf = 1

    if (grepl("wide", subtype)) {
      sampler$tau2$sc_adj <- 1.5
    } else {
      sampler$tau2$sc_adj <- 1.0
    }

  }

}

(state <- mc_out$state)

### timing run
mc_time <- time_hs(state = state, prior = prior, data = dat,
                   sampler = sampler, n_iter = n_iter, param = "tau2",
                   upper_tau2 = 1.0e9)
                 # , appx_threshold = 0.0)

mc_time$timing
(samp_p_sec <- mc_time$timing$EffSamp / mc_time$timing$userTime)


### diagnostics
draws_tau2 <- mc_time$draws
(tau2_SE <- summary(as.mcmc(draws_tau2))$statistics["Time-series SE"])
(tau2_mn <- mean(draws_tau2))
(tau2_sd <- sd(draws_tau2))

# plot(as.mcmc(draws_tau2))
# plot(as.mcmc(log(draws_tau2)))

if (type == "Qslice") {
  mc_out <- mcmc_hs(state = state, prior = prior, data = dat,
                        sampler = sampler,
                        n_iter = 2e3,
                        save = TRUE, prog = 0,
                        upper_tau2 = 1.0e9)
                        # , appx_threshold = 0.0) # samples for pseudo-target

  draws_u <- sapply(mc_out$extras, function(x) x$u)
  (AUC <- auc(u = draws_u))
} else {
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
                     nEval = mc_time$timing$nEval,
                     ESS = mc_time$timing$EffSamp,
                     userTime = mc_time$timing$userTime,
                     sysTime = mc_time$timing$sysTime,
                     elapsedTime = mc_time$timing$elapsedTime,
                     sampPsec = samp_p_sec,
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
