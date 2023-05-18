## A script to summarize the gprior results
# Author: Sam Johnson

library(magrittr)
source('formatingFunctions.R')

# reading in all the files
files <- (Sys.glob("data/*.rds"))

data <- sapply(files, FUN = \(files) readRDS(files), simplify = FALSE)

results <- lapply(data, resultsFunc)

gResults <- sapply(1:length(results), FUN = \(i) {
  name <- names(results)[i]
  name <- stringr::str_remove(name, pattern = 'data/')
  name <- stringr::str_remove(name, pattern = '.rds')
  list <- results[[i]]
  tempDf <- list[1,]
  tempDf$method <- name
  
  tempDf
})


t(gResults) %>% 
  data.frame() %>% 
  dplyr::relocate(method, .before = 'lwrCrd')


