# setwd('~/cucumber/sam_comparison/curve1')
source('norm_setup.R')

##
#### Stepping out ####
##

load('input/steppingout.rds')

i <- commandArgs(trailingOnly = TRUE)[1]

row <- trials_stepping_out[i,]
metrics <- stepping_out_time_eval(
  samples = row$samples[1],
  x_0 = row$x[1],
  lf_func = lf,
  w_value = row$w[1],
  max_value = Inf,
  log_value = TRUE
)
thin <- (length(metrics$Draws[[1]])/metrics$EffSamp) * 10
thinDraws <- LaplacesDemon::Thin(metrics$Draws[[1]], thin)
tempDf <- data.frame(nEval = metrics$nEval,
           ESS = metrics$EffSamp,
           userTime = metrics$userTime,
           sysTime = metrics$sysTime,
           elapsedTime = metrics$elapsedTime,
           sampPsec = metrics$EffSamp/metrics$userTime,
           w = row$w,
           ksTest = ks.test(thinDraws, pnorm, 0, 1)$p.value)

# saving samples
write.table(tempDf, file = paste0('output/steppingout/steppinout',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')