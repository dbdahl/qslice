rm(list=ls())
library("tidyverse")

set.seed(240329)

n_reps <- 100

targets <- c("hyper-g", "hyper-g-log")

types <- c("rw", "stepping", "latent", "gess", "imh", "Qslice")
subtypes <- c("AUC_samples", "Laplace_analytic", "Laplace_analytic_wide")

samplers0 <- rbind( data.frame(type = c("rw", "stepping", "latent"), subtype = NA),
                   expand.grid(subtype = subtypes, type = c("gess", "imh", "Qslice"))[,2:1]
)
samplers0

samplers1 <- lapply(1:length(targets), function(i) cbind(target = targets[i], samplers0)) %>% do.call(rbind, .)
samplers1

sched <- lapply(1:n_reps, function(i) cbind(samplers1, rep = i)) %>% do.call(rbind, .)
sched

(n_jobs <- nrow(sched))

job_order <- sample(n_jobs, size = n_jobs, replace = FALSE)

save(file = paste0("schedule_all.rda"),
     sched, n_jobs, job_order)
