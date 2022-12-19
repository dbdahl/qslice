# setwd('~/cucumber/sam_comparison/curve2')
source('curve2_setup.R')

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
  dplyr::mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = length(draws)/ESS,
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    truncThinDraws = list(
      sample(unlist(thinDraws), ifelse(samplesThin <= sampleSize, samplesThin, sampleSize))
    ),
    ksTest = ks.test(truncThinDraws, cdf)$p.value
  ) %>%
  dplyr::select(-c(metrics, truncThinDraws)) %>%
  dplyr::mutate(SampPSec = ESS / time) %>%
  dplyr::relocate(samples, .after = time)

# # the functions to compare against
# dist_df <- data.frame(
#   dist = rep("cdf", 6),
#   mean2 = c(20, 20, 20, 19, 21, 20),
#   sd2 = c(4, 2, 3, 1, 1, 1)
# )
# 
# list_hldr <- list(length = nrow(dist_df))
# for (i in 1:nrow(dist_df)) {
#   list_hldr[[i]] <- lapply(
#     rand_walk_metrics$thinDraws,
#     FUN = ks.test,
#     y = cdf,
#     mean2 = dist_df[i, 'mean2'],
#     sd2 = dist_df[i, 'sd2']
#   )
# }
# 
# 
# # saving the power test
# pdf(file = "../../images_slice_sampler_comp/curve2_rand_walk_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(rand_walk_metrics,paste0("../../data/curve",curve_num,"_rand_walk_metrics"))

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "../../images_slice_sampler_comp/curve2_rand_walk.pdf")
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2
)
lapply(rand_walk_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black'))
})
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2,
  add = TRUE
)
dev.off()


rand_walk_min_max <-
  results(rand_walk_metrics, method = "Random Walk")

