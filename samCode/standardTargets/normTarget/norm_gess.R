source('norm_setup.R')

##
#### gess ####
##

load('input/gess.rds')

i <- commandArgs(trailingOnly = TRUE)[1]

row <- trials_gess[i,]
metrics <- gess_time_eval(
  samples = row$samples[1],
  x_0 = row$x[1],
  lf_func = lf,
  mu_value = row$mu[1],
  sigma_value = row$sigma[1],
  df_value = row$df[1],
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
           mu = row$mu,
           sigma = row$sigma,
           df = row$df,
           ksTest = ks.test(thinDraws, pnorm, 0, 1)$p.value)

# saving samples
write.table(tempDf, file = paste0('output/gess/gess',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
