# how many lines are in a rds file
# authors: Sam Johnson, Matt Heiner

args <- commandArgs(trailingOnly = TRUE)
data_use <- args[1]

file <- paste0("schedule_all_", data_use, ".rda")

load(file)

cat(n_jobs)

quit(save = "no")
