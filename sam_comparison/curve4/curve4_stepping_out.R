# setwd('~/cucumber/sam_comparison/curve4')
source('curve4_setup.R')

##
#### Stepping out ####
##

# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(
    metrics = stepping_out_time_eval(
      samples = samples,
      x_0 = x,
      lf_func = lf,
      w_value = w,
      max_value = Inf,
      log_value = TRUE
    )
  ) %>%
  mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = (length(draws)/ESS) * 10,
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, cdf)$p.value
  ) %>%
  dplyr::select(-metrics, -thinDraws) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
saveRDS(stepping_out_metrics,paste0("../../data/curve",curve_num,"_stepping_out_metrics"))

