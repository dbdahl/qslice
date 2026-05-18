rm(list=ls())
library("tidyverse")

set.seed(260518)

n_reps <- 50

data_use <- "db"
# data_use <- "db40"

targets <- paste(data_use, c("tau2_marg", "tau2_marg-log"), sep = "_")

types <- c("rw", "stepping", "latent", "gess", "imh", "Qslice")
subtypes <- c("AUC_samples", "MSW_samples") #, "Laplace_analytic", "Laplace_analytic_wide") # Add "MSW_samples"?

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

save(file = paste0("schedule_all_", data_use, ".rda"),
     sched, n_jobs, job_order)

