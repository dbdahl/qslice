# how many lines are in a rds file
# author: Sam Johnson

file <- commandArgs(trailingOnly = TRUE)

load(file)

dataName <- ls()[2]

numOfLines <- if (dataName == 'trials_gess') {
  nrow(trials_gess)
} else if (dataName == 'trials_stepping_out') {
  nrow(trials_stepping_out)
} else if (dataName == 'trials_transform') {
  nrow(trials_transform)
} else if (dataName == 'trials_rand_walk') {
  nrow(trials_rand_walk)
} else if (dataName == 'trials_latent') {
  nrow(trials_latent)
}

cat(numOfLines)