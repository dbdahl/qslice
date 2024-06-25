args <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(args[1])  # job id
dte <- as.numeric(args[2])

##### for testing
# ii <- 1944
# dte <- 240520
#####

library("qslice")
library("coda")
source("0_data.R")
source("0_prior.R")
source("functions/MH_samplers.R")
source("functions/full_conditionals.R")
source("functions/mcmc_gprior.R")
source("functions/tune.R")
source("functions/Laplace_g.R")
source("functions/invgamma.R")

load("schedule_all.rda")

n_iter <- 50e3

run_id <- job_order[ii]
(run_info <- sched[run_id,])
target <- as.character(run_info$target)
logG <- grepl("log", target)
type <- as.character(run_info$type)
subtype <- as.character(run_info$subtype)

state <- list(iter = 0)
state$psi <- prior$a_0 / prior$b_0
state$g <- prior$g_max / 3
state$beta <- prior$beta_0
state

sampler <- list()

sampler$g <- list(type = "stepping", subtype = NA, logG = logG,
                  w = ifelse(logG, 5.0, 50.0))


### initial burn-in
mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = 8e3,
                      save = FALSE, prog = 0)

(state <- mc_out$state)

mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                      sampler = sampler,
                      n_iter = 2e3,
                      save = TRUE, prog = 0) # samples for pseudo-target

(state <- mc_out$state)


### Additional prep
if (type %in% c("rw", "stepping", "latent")) { # will require tuning

  if (type == "latent") {
    tune_bnds_init = c(0.001, 0.5)
    state$latent_s <- ifelse(logG, 10.0, 100.0)
  } else {
    if (isTRUE(logG)) {
      tune_bnds_init = c(0.5, 5.0)
    }
    tune_bnds_init = c(5.0, 100.0)
  }

  sampler$g$type <- type
  sampler$g$subtype <- subtype

  mc_tune <- tune(state = state, prior = prior, data = dat,
                  sampler = sampler,
                  bnds_init = tune_bnds_init,
                  n_iter = 1e3, n_grid = 5, n_rep = 3,
                  n_rounds = 5, range_frac = 0.5,
                  verbose = TRUE)

  sampler$g <- mc_tune$sampler$g

  if (type == "latent") {
    state <- mc_tune$state # since latent_s was burned in during tuning
  }

} else if (grepl("samples", subtype)) { # tune with samples

  mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                        sampler = sampler,
                        n_iter = 1e3,
                        save = TRUE, prog = 0)

  g_samples <- sapply(mc_out$sims, function(x) x$g)
  if (isTRUE(logG)) {
    g_samples <- log(g_samples)
  }

  tmp_pseu <- pseudo_opt(samples = g_samples,
                         type = "samples",
                         family = "t",
                         degf = c(1, 5),
                         lb = ifelse(logG, -Inf, 0.0),
                         ub = ifelse(logG, log(prior$g_max), prior$g_max),
                         utility_type = "AUC",
                         plot = FALSE)

  sampler$g$type <- type
  sampler$g$subtype <- subtype
  sampler$g$pseudo <- tmp_pseu$pseudo
  sampler$g$loc <- tmp_pseu$pseudo$params$loc
  sampler$g$sc <- tmp_pseu$pseudo$params$sc
  sampler$g$degf <- tmp_pseu$pseudo$params$degf
  sampler$g$txt <- tmp_pseu$pseudo$txt

} else if (grepl("Laplace_analytic", subtype)) {

  sampler$g$type <- type
  sampler$g$subtype <- subtype

  if (isTRUE(logG)) {

    sampler$g$degf = 5

    if (grepl("wide", subtype)) {
      sampler$g$sc_adj <- 1.2
    } else {
      sampler$g$sc_adj <- 1.0
    }

  } else {

    sampler$g$degf = 1

    if (grepl("wide", subtype)) {
      sampler$g$sc_adj <- 1.5
    } else {
      sampler$g$sc_adj <- 1.0
    }

  }

}

### timing run
mc_time <- time_gprior(state = state, prior = prior, data = dat,
                       sampler = sampler, n_iter = n_iter)

mc_time$timing
(samp_p_sec <- mc_time$timing$EffSamp / mc_time$timing$userTime)


### diagnostics
draws_g <- mc_time$draws
(g_SE <- summary(as.mcmc(draws_g))$statistics["Time-series SE"])
(g_mn <- mean(draws_g))
(g_sd <- sd(draws_g))

if (type == "Qslice") {
  mc_out <- mcmc_gprior(state = state, prior = prior, data = dat,
                        sampler = sampler,
                        n_iter = 2e3,
                        save = TRUE, prog = 0) # samples for pseudo-target

  draws_u <- sapply(mc_out$extras, function(x) x$u)
  # hist(draws_u, breaks = 30)
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
                     g_mn = g_mn,
                     g_sd = g_sd,
                     g_SE = g_SE,
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
