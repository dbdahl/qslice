# setwd('~/cucumber/sam_comparison/curve7')
source('curve7_setup.R')

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
    thin = length(draws)/ESS,
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    truncThinDraws = list(
      sample(unlist(thinDraws), ifelse(samplesThin <= sampleSize, samplesThin, sampleSize))
    ),
    ksTest = ks.test(truncThinDraws, pnorm, mean = 20, sd = 5)$p.value
  ) %>%
  select(-metrics, -truncThinDraws) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# # the functions to compare against
# dist_df <- data.frame(
#   dist = rep("pnorm", 6),
#   mean = c(18, 21, 19, 20, 20, 20),
#   sd = c(5, 5, 5, 6, 4, 5)
# )
# 
# list_hldr <- list(length = nrow(dist_df))
# for (i in 1:nrow(dist_df)) {
#   list_hldr[[i]] <- lapply(
#     stepping_out_metrics$thinDraws,
#     FUN = ks.test,
#     y = dist_df[i, 'dist'],
#     mean = dist_df[i, 'mean'],
#     sd = dist_df[i, 'sd']
#   )
# }
# 
# # saving the power test
# pdf(file = "../../images_slice_sampler_comp/curve7_stepping_out_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out the metrics table
saveRDS(stepping_out_metrics,paste0("../../data/curve",curve_num,"_stepping_out_metrics"))

pdf(file = "../../images_slice_sampler_comp/curve7_stepping_out.pdf")
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2
)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
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

stepping_min_max <-
  results(stepping_out_metrics, method = "Stepping Out")

rm(
  # list_hldr,
  # dist_df,
  trials_stepping_out,
  stepping_out_metrics
)
