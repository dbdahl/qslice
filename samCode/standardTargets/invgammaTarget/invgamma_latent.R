# setwd('~/cucumber/sam_comparison/curve1')
source('invgamma_setup.R')

##
#### latent ####
##

load('input/latent.rds')

i <- commandArgs(trailingOnly = TRUE)[1]

row <- trials_latent[i,]
metrics <- latent_time_eval(
  samples = row$samples[1],
  x_0 = row$x[1],
  lf_func = lf,
  s_0 = row$s[1],
  rate_value = row$rate[1],
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
                     s = row$s,
                     rate = row$rate,
                     ksTest = ks.test(thinDraws, pinvgamma, 2, 1)$p.value)

# saving samples
write.table(tempDf, file = paste0('output/latent/latent',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')