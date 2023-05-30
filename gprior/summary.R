## A script to summarize the gprior results
# Author: Sam Johnson

library(magrittr)
source('formatingFunctions.R')

# reading in all the files
files <- (Sys.glob("data/*.rds"))

data <- sapply(files, FUN = \(files) readRDS(files), simplify = FALSE)

results <- lapply(data, resultsFunc)

# getting the g results

gResults <- sapply(1:length(results), FUN = \(i) {
  name <- names(results)[i]
  name <- stringr::str_remove(name, pattern = 'data/')
  name <- stringr::str_remove(name, pattern = '.rds')
  list <- results[[i]]
  tempDf <- list[1,]
  tempDf$method <- name
  
  tempDf
})


temp <- t(gResults) %>% 
  data.frame() %>% 
  dplyr::relocate(method, .before = 'lwrCrd') %>% 
  dplyr::arrange(desc(unlist(SampPSec)))

temp

## plotting the densities

gDensities <- lapply(data, FUN = \(list) {
  g <- sapply(1:length(list), FUN = \(i) {
    list[[i]]$g
  })
  c(g)
})

colors <- RColorBrewer::brewer.pal(length(gDensities), name = 'Paired')
names <- names(gDensities) %>% stringr::str_remove('data/') %>% stringr::str_remove('.rds')

plot(density(gDensities[[1]]), col = colors[1], main = 'Density Plots of g', ylim = c(0,0.05), lty = 1)
for(i in 2:length(gDensities)) {
  lines(density(gDensities[[i]]), col = colors[i], lty = i)
}
legend(x = 60, y = 0.048, legend = names, col = colors, lty = 1:length(gDensities), cex = 0.8)


## getting the info on the u's
transformData <- sapply(files[grepl('*transform*',files)],  FUN = \(files) readRDS(files), simplify = FALSE)

sapply(1:length(transformData), FUN = \(i) {
  list <- transformData[[i]]
  u <- c(sapply(list, FUN = \(x) x$u))
  aucDiagnostic <- auc_diagnostic(u)
  hist(u, main = names(transformData)[i], sub = paste0('auc diag:', round(aucDiagnostic)))
})


apply(sapply(transformData[[1]], FUN = \(list) list$u),2,hist)
