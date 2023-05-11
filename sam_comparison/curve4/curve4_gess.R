# setwd('~/cucumber/sam_comparison/curve4')
source('curve4_setup.R')

##
#### gess ####
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename(
    'samples' = 'Var1',
    'x' = 'Var2',
    'mu' = 'Var3',
    'sigma' = 'Var4',
    'df' = 'Var5'
  )

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(
    metrics = gess_time_eval(
      samples = samples,
      x_0 = x,
      mu_value = mu,
      sigma_value = sigma,
      df_value = df,
      lf_func = lf,
      lwr_support = 0,
      log_value = TRUE
    )
  )  %>%
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
  select(-metrics, -thinDraws) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
saveRDS(gess_metrics,paste0("../../data/curve",curve_num,"_gess_metrics"))

