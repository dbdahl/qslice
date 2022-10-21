source('curve7_setup.R')

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
    thin = min(which(
      acf(draws, plot = FALSE, lag.max = 5000)$acf < auto.cor.lim
    )),
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, pnorm, mean = 20, sd = 5)$p.value
  ) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# the functions to compare against
dist_df <- data.frame(
  dist = rep("pnorm", 6),
  mean = c(18, 21, 19, 20, 20, 20),
  sd = c(5, 5, 5, 6, 4, 5)
)

list_hldr <- list(length = nrow(dist_df))
for (i in 1:nrow(dist_df)) {
  list_hldr[[i]] <- lapply(
    latent_metrics$thinDraws,
    FUN = ks.test,
    y = dist_df[i, 'dist'],
    mean = dist_df[i, 'mean'],
    sd = dist_df[i, 'sd']
  )
}

# saving the power test
pdf(file = "../../images_slice_sampler_comp/curve7_latent_power.pdf")
extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
dev.off()

# printing out the metrics table
saveRDS(latent_metrics,paste0("../../data/curve",curve_num,"_latent_metrics"))

# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "../../images_slice_sampler_comp/curve7_latent.pdf")
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2
)
lapply(latent_metrics$thinDraws, function(x) {
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

latent_min_max <- results(latent_metrics, method = "Latent")

rm(
  list_hldr,
  dist_df,
  trials_latent,
  latent_metrics
)