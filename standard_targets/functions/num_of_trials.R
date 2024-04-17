# how many lines are in a rds file
# author: Sam Johnson

args <- commandArgs(trailingOnly = TRUE)

target <- args[1]
rnd <- args[2] |> as.numeric()

file <- paste0("input/schedule_target", target, "_round", rnd, ".rda")

load(file)

cat(n_jobs)

quit(save = "no")
