args <- commandArgs(trailingOnly = TRUE)
rnd <- as.numeric(args[1])  # round (1 for tuning, 2 for final)
ii <- as.numeric(args[2])  # job id
dte <- as.numeric(args[3])

##### for testing
# rnd <- 2
# ii <- 1935
# dte <- 240520
#####

source("0_setup.R")

load(paste0("input/schedule_all_round", rnd, ".rda"))

run_id <- job_order[ii]
(run_info <- sched[run_id,])
target <- as.character(run_info$target)
type <- as.character(run_info$type)
subtype <- as.character(run_info$subtype)

source(paste0("0_setup_", target, ".R"))
load(paste0("input/schedule_target", target, "_round", rnd, ".rda")) # overwrites previous sched, so don't use sched or job order again

n_runs <- 2 # more than 2 if running on Intel machine with MKL

if (isTRUE(run_info["subtype"] %in% c("MSW_samples", "AUC_samples", "MM_Cauchy"))) {
  set.seed(dte + run_info$rep)
  source("1_pseudo_from_samples.R")  # pseudo-target from samples
  run_info$algo_descrip = trials[[type]][[subtype]]$algo_descrip
}

set.seed(dte + ii)

if (type %in% c("Qslice", "imh")) {
  settings <- trials[[type]][[subtype]]$pseudo$pseu
  n_iter <- trials[[type]][[subtype]]$n_iter
  x0 <- trials[[type]][[subtype]]$x0
} else {
  settings <- trials[[type]]$params[run_info$setting, , drop = FALSE]
  n_iter <- trials[[type]]$n_iter
  x0 <- trials[[type]]$x0
}

for (jj in 1:n_runs) {
  mcmc_time_out <- sampler_time_eval(
    type = type,
    n_iter = n_iter,
    lf_func = truth$ld,
    support = c(truth$lb, truth$ub),
    x_0 = x0,
    settings = settings
  )
}

thin <- min(10 * n_iter / mcmc_time_out$EffSamp, 100)
thinDraws <- LaplacesDemon::Thin(mcmc_time_out$Draws[[1]], thin)

tempDf <- data.frame(round = rnd,
                     target = target,
                     type = type,
                     subtype = subtype,
                     setting = run_info$setting,
                     rep = run_info$rep,
                     run_id = run_id,
                     ii = ii,
                     algo_descrip = as.character(run_info$algo_descrip),
                     n_iter = n_iter,
                     x_0 = x0,
                     nEval = mcmc_time_out$nEval,
                     ESS = mcmc_time_out$EffSamp,
                     userTime = mcmc_time_out$userTime,
                     sysTime = mcmc_time_out$sysTime,
                     elapsedTime = mcmc_time_out$elapsedTime,
                     sampPsec = mcmc_time_out$EffSamp / mcmc_time_out$userTime,
                     ks_pval = ks.test(thinDraws, truth$p)$p.value)

# saving samples
write.table(tempDf, file = paste0("output/round", rnd,
                                  "_target", target,
                                  "_type", type,
                                  "_subtype", subtype,
                                  "_setting", run_info$setting,
                                  "_rep", run_info$rep,
                                  "_dte", dte,
                                  ".csv"),
            append = FALSE, sep = ",",
            row.names = FALSE, col.names = TRUE)

print('finished')

quit(save = "no")
