rm(list=ls())

set.seed(230308)

targets <- c("normal", "gamma", "igamma")
targets <- c("normal", "gamma", "igamma", "gammalog", "igammalog")

rnd <- 2

sched_list <- list()

for (tg in targets) { # requires 1_setup_trials.R be run for each target first
  load(paste0("input/schedule_target", tg, "_round", rnd, ".rda"))
  sched_list[[tg]] <- sched
}
rm(sched)

sched <- do.call(rbind, sched_list)

rm(list = setdiff(ls(), c("targets", "rnd", "sched")))
ls()

str(sched)
head(sched); tail(sched)

(n_jobs <- nrow(sched))

job_order <- sample(n_jobs, size = n_jobs, replace = FALSE)

save(file = paste0("input/schedule_all_round", rnd, ".rda"),
     rnd, sched, n_jobs, job_order)
