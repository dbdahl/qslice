# setwd('~/cucumber/sam_comparison/curve1')
source('norm_setup.R')

##
#### latent ####
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename(
    'samples' = 'Var1',
    'x' = 'Var2',
    's' = 'Var3',
    'rate' = 'Var4'
  )

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(
    metrics = latent_time_eval(
      samples = samples,
      x_0 = x,
      s_0 = s,
      lf_func = lf,
      rate_value = rate,
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
    ksTest = ks.test(thinDraws, pnorm, 0, 1)$p.value
  ) %>%
  select(-metrics, -thinDraws) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

if(saveInd) saveRDS(latent_metrics, file = 'data/latent.rds')

