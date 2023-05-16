## A script to summarize the gprior results
# Author: Sam Johnson

source('formatingFunctions.R')

# reading in all the files
files <- (Sys.glob("data/*.rds"))

data <- sapply(files, FUN = \(files) readRDS(files), simplify = FALSE)

results <- lapply(data, resultsFunc)

