source('curve6_setup.R')

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
      upr_support = 1,
      lwr_support = 0,
      log_value = TRUE
    )
  )  %>%
  mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = min(which(
      acf(draws, plot = FALSE, lag.max = 10000)$acf < auto.cor.lim
    )),
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, pbeta, shape1 = 0.2, shape2 = 0.8)$p.value
  ) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# the functions to compare against
dist_df <- data.frame(
  dist = rep("pbeta", 6),
  shape1 = c(.1, .8, .5, .2, .2, .2),
  shape2 = c(.9, .8, .8, .5, .2, .8)
)

list_hldr <- list(length = nrow(dist_df))
for (i in 1:nrow(dist_df)) {
  list_hldr[[i]] <- lapply(
    gess_metrics$thinDraws,
    FUN = ks.test,
    y = dist_df[i, 'dist'],
    shape1 = dist_df[i, 'shape1'],
    shape2 = dist_df[i, 'shape2']
  )
}

# saving the power test
pdf(file = "../../images_slice_sampler_comp/curve6_gess_power.pdf")
extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
dev.off()

# printing out the metrics table
saveRDS(gess_metrics,paste0("../../data/curve",curve_num,"_gess_metrics"))

pdf(file = "../../images_slice_sampler_comp/curve6_gess.pdf")
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2
)
lapply(gess_metrics$thinDraws, function(x) {
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


gess_min_max <- results(gess_metrics, method = "GESS")

rm(
  list_hldr,
  dist_df,
  trials_gess,
  gess_metrics
)
