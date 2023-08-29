## simulation study on the normal curve
## author: Sam Johnson

source('norm_setup.R')

##
#### Transform ####
##

load('input/transform.rds')

i <- commandArgs(trailingOnly = TRUE)[1]

row <- trials_transform[i,]
metrics <- transform_time_eval(
  samples = row$samples[1],
  x_0 = row$x[1],
  lf_func = lf,
  pseudo_pdf_log = row$log_pdf[[1]],
  pseudo_cdf_inv = row$inv_cdf[[1]],
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
           t = row$t[[1]],
           ksTest = ks.test(thinDraws, pnorm, 0, 1)$p.value)

# saving samples
write.table(tempDf, file = paste0('output/transform/transform',i,'.csv'), append = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

print('finished')
