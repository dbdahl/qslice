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


temp <- t(gResults) %>% 
  data.frame() %>% 
  dplyr::relocate(method, .before = 'lwrCrd') %>% 
  dplyr::arrange(desc(unlist(SampPSec)))

temp

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
