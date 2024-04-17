# how many lines are in a rds file
# author: Sam Johnson

args <- commandArgs(trailingOnly = TRUE)

rnd <- args[1] |> as.numeric()

file <- paste0("input/schedule_all_round", rnd, ".rda")

load(file)

cat(n_jobs)

quit(save = "no")
