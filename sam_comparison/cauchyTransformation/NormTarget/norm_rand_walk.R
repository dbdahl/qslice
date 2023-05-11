# setwd('~/cucumber/sam_comparison/curve1')
source('norm_setup.R')

##
#### Metropolis Random Walk ####
##

trials_rand_walk <- expand.grid(samples, x, c) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'c' = 'Var3')

# evaluating the sampler at each different iteration
rand_walk_metrics <- trials_rand_walk %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    metrics = random_walk_time_eval(samples = samples,
                                    x_0 = x,
                                    c = c,
                                    lf_func = lf)
  ) %>%
  mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    acceptanceRate = metrics$acceptance_rate,
    thin = (length(draws)/ESS) * 10,
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, pnorm, 0, 1)$p.value
  ) %>%
  dplyr::select(-metrics, -thinDraws) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

if(saveInd) saveRDS(rand_walk_metrics, file = 'data/randWalk.rds')