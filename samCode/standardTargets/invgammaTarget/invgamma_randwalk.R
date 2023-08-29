# setwd('~/cucumber/sam_comparison/curve1')
source('invgamma_setup.R')

##
#### Metropolis Random Walk ####
##

load('input/randwalk.rds')

i <- commandArgs(trailingOnly = TRUE)[1]

row <- trials_rand_walk[i,]
metrics <- random_walk_time_eval(
  samples = row$samples[1],
  x_0 = row$x[1],
  lf_func = lf,
  c = row$c[1]
)
thin <- (length(metrics$Draws[[1]])/metrics$EffSamp) * 10
thinDraws <- LaplacesDemon::Thin(metrics$Draws[[1]], thin)
tempDf <- data.frame(nEval = metrics$nEval,
                     ESS = metrics$EffSamp,
                     userTime = metrics$userTime,
                     sysTime = metrics$sysTime,
                     elapsedTime = metrics$elapsedTime,
                     sampPsec = metrics$EffSamp/metrics$userTime,
                     c = row$c,
                     acceptanceRate = metrics$acceptanceRate,
                     ksTest = ks.test(thinDraws, pinvgamma, 2, 1)$p.value)

# saving samples
write.table(tempDf, file = paste0('output/randwalk/randwalk',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
