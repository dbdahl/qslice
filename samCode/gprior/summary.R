## A script to summarize the gprior results
# Author: Sam Johnson

# sourcing the files the different samplers
filesToSource <- Sys.glob('../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)
library(magrittr)
library(tidyverse)
library(forcats)
source('formatingFunctions.R')
source('setupGprior.R')

theme_set(
  theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = '',
          plot.title = element_text(hjust = 0.5, size = 30),
          axis.title.x = element_text(size = 25), 
          axis.text = element_text(size = 20),
          panel.grid.major = element_blank()
          # axis.line.y.left = element_line(color = 'black'),
          # axis.line.x.bottom = element_line(color = 'black'),
          # axis.ticks = element_line(color = 'black'),
          # axis.ticks.length=unit(.25, "cm"),
          # plot.title = element_text(hjust = 0.5)
          
    )
)

# reading in all the files
# files <- (Sys.glob("data/*.rds"))
# data <- sapply(files, readRDS, simplify = FALSE)
# results <- lapply(data, resultsFunc)

files <- Sys.glob("output/*/*.csv")
data <- sapply(files, read.csv, simplify = FALSE)
results <- do.call(rbind, data)  |> 
  mutate(
    sampler = case_when(
      grepl('*Samples*|*auto*|*laplace*',t) ~ 'Transform',
      grepl('*latent*',t) ~ 'Latent',
      grepl('*gess*',t) ~ 'GESS',
      grepl('*steppingout*',t) ~ 'Stepping Out',
      grepl('*independence*',t) ~ 'Independence',
      grepl('*randwalk*',t) ~ 'Random Walk',
      TRUE ~ 'Other'
    ),
    pseudoType = case_when(
      grepl('*optimSamplesAUC*',t) ~ 'OSAUC',
      grepl('*optimSample*',t) ~ 'OS',
      grepl('*auto*',t) ~ 'Auto',
      grepl('*laplace*',t) ~ 'Laplace',
      TRUE ~ 'Other'
    ))

pdf(file = 'images/boxplot.pdf')
results |> 
  ggplot(aes(x = sampPsec, y = sampler, fill = t)) +
  geom_boxplot() +
  theme(
    legend.position = 'bottom'
  ) +
  labs(
    title = 'Zellner G Prior',
    y = '',
    x = 'Effective Samples Per CPU Second',
    fill = 'Pseudo Target Method'
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.title.x = element_text(size = 20)
  )
dev.off()

pdf(file = 'images/bestboxplot.pdf')
results |> 
  dplyr::filter(t %in% c('optimSamplesAUC','latent','steppingout','gess','independence','randwalk')) |> 
  ggplot(aes(x = sampPsec, y = sampler, fill = sampler)) +
  geom_boxplot() +
  labs(
    title = 'Zellner G Prior',
    y = '',
    x = 'Effective Samples Per CPU Second'
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.title.x = element_text(size = 20)
  )
dev.off()


# transform
gSamples <- gprior_sampler(list(Nsamples = 5000, Nburnin = 100,Nthin = 1, method = 'SteppingOut', w = 15))
# optim samples
optimSamplesAUCList$approxSamples <- gSamples$g
# chainSamplesOptimSamples <- sapply(1:nChains, FUN = \(i) { gprior_sampler(optimSamplesAUCList) } )
samplesOptimSamplesAUC <- gprior_sampler(optimSamplesAUCList)

pdf('images/uhist.pdf')
data.frame(u = samplesOptimSamplesAUC$u) |> 
  ggplot(aes(x = u)) +
  geom_histogram(bins = 30) +
  labs(
    title = paste0('Water: ', round(water_area(u = samplesOptimSamplesAUC$u),2))
  )
dev.off()

auc(u = samplesOptimSamplesAUC$u)
water_area(u = samplesOptimSamplesAUC$u)

# function to find how much water a curve can hold
water_area <- function(x, y, u = NULL, nbins = 30, plot = FALSE) {
  # browser()
  if(is.null(u)) {
    
    if (is.function(y)) {
      y <- y(x)
    } 
    
    stopifnot (length(y) == length(x))
    stopifnot(all(x > 0.0) && all(y >= 0.0) && all(x < 1.0))
    
  } else {
    
    bins <- seq(0.0, 1.0, len = nbins + 1)
    y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins=nbins)
    x <- (bins[-(nbins+1)] + bins[-1]) / 2.0
    
  }
  
  n <- length(y)
  
  # getting a normalized y
  yn <- y / max(y)
  # creating a vector to store which indecies are underwater
  underwater <- vector(length = n)
  # finding which values are underwater
  for( i in 2:(n-1) ) {
    if( max(yn[1:(i-1)]) > yn[i] && max(yn[(i+1):n]) > yn[i] ) {
      underwater[i] <- TRUE
    }
  }
  # checking if any points are underwater
  if( sum(underwater) == 0 ) {
    if( plot ) plot(x,yn,type = 'l')
    return(0)
  }
  # finding which indicies are under water
  underwaterIndex <- (1:n)[underwater]
  # finding the different bowls
  splits <- which(diff(underwaterIndex) != 1) + 1
  bowls <- unname(split(underwaterIndex, cumsum(seq_along(underwaterIndex) %in% splits)))
  # finding the areas of the bowls
  bowlArea <- sapply(1:length(bowls), FUN = \(i) {
    list <- bowls[[i]]
    waterLevel <- min(yn[c((min(list)-1),(max(list)+1))])
    area <- (waterLevel - yn[list]) |> sum()
    c(area = area, waterLevel = waterLevel)
  })
  # plotting the curve
  if( plot ) {
    plot(x,yn,bty = 'l', type = 'l', main = '', ylab = '', xlab = '')
    # points(x[underwater], yn[underwater], col = 'red')
    lapply(bowls, FUN = \(list) {
      for( i in 1:length(list) ) {
        index <- list[i]
        x0 <- x[index]
        y0 <- yn[index]
        y1 <- min(yn[c((min(list)-1),(max(list)+1))])
        segments(x0 = x0, x1 = x0, y0 = y0, y1 = y1, col = 'dodgerblue2')
      }
    })
  }
  sum(bowlArea[1,]) / n
}

pdf(file = 'images/sidebyside.pdf')
par(mfrow = c(1,2))
hist(samplesOptimSamplesAUC$u, main = '', xlab = '', bty = 'l',)
water_area(u = samplesOptimSamplesAUC$u, plot = TRUE)
dev.off()