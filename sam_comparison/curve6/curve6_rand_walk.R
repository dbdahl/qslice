# setwd('~/cucumber/sam_comparison/curve6')
source('curve6_setup.R')

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
                                    lf_func = lf,
                                    support = c(0.0001, .9999))
  ) %>%
  mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = (length(draws)/ESS) * 10,
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws,  pbeta, shape1 = 0.2, shape2 = 0.8)$p.value
  ) %>%
  dplyr::select(-c(metrics,thinDraws)) %>%
  dplyr::mutate(SampPSec = ESS / time) %>%
  dplyr::relocate(samples, .after = time)

# printing out metrics table
saveRDS(rand_walk_metrics,paste0("../../data/curve",curve_num,"_rand_walk_metrics"))

