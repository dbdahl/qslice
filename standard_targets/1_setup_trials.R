rm(list=ls())

#### inputs

target <- "normal"
target <- "gamma"
target <- "gammalog"
target <- "igamma"
target <- "igammalog"

rnd <- 1 # round (1 for exploring settings; 2 for final comparisons at selected settings)
n_rep <- 5 # number of replicate runs for each sampler/setting

x0 <- 0.2  # same initial value for all samplers
wide_factor <- 4.0 # scale inflation for methods using a "diffuse" pseudo-target

#########

library("qslice")

source(paste0("0_setup_", target, ".R"))
source("1_algo_params.R")

n_iter <- 50e3
trials <- list()

competitors <- c("rw", "stepping", "gess", "latent")

for (cc in competitors) {
  trials[[cc]] <- list(target = target,
                       type = cc,
                       subtype = NA,
                       n_iter = n_iter,
                       x0 = x0,
                       params = algo_params[[target]][[cc]])
}


### create schedule
sched0 <- list()


if (rnd > 1) {
  source("1_setup_trials_pseudo.R") # methods with pseudo-targets won't require initial calibration
}


for (cc in competitors) {

  sched0[[cc]] <- expand.grid(target = target,
                              type = cc, subtype = NA,
                              setting = 1:nrow(trials[[cc]]$params),
                              rep = 1:n_rep,
                              stringsAsFactors = FALSE)
  sched0[[cc]]$algo_descrip <- sapply(1:nrow(sched0[[cc]]), function(i) {
    paste(colnames(trials[[cc]]$params),
          trials[[cc]]$params[sched0[[cc]]$setting[i],], collapse = "; ")
  })

}

sched <- do.call(rbind, sched0)
rownames(sched) <- NULL
head(sched); tail(sched)

(n_jobs <- nrow(sched))
job_order <- sample(n_jobs, size = n_jobs, replace = FALSE)

sched[job_order,] |> head(n = 30)
sched[job_order,] |> tail(n = 30)

save(file = paste0("input/schedule_target", target, "_round", rnd, ".rda"),
     target, rnd, n_rep, x0, n_iter, trials, sched, n_jobs, job_order)

