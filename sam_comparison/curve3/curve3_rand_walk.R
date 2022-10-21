source('curve3_setup.R')

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
                                    support = c(0, Inf))
  ) %>%
  dplyr::mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = min(which(
      acf(draws, plot = FALSE, lag.max = 5000)$acf < auto.cor.lim
    )),
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, pgamma, shape = 2.5, rate = 1)$p.value
  ) %>%
  dplyr::select(-c(metrics)) %>%
  dplyr::mutate(SampPSec = ESS / time) %>%
  dplyr::relocate(samples, .after = time)

# the functions to compare against
dist_df <- data.frame(
  dist = rep("pgamma", 6),
  shape = c(1, 2, 2.5, 5, 3, 2.5),
  rate = c(1, 1, 5, 1, 1, 1)
)

list_hldr <- list(length = nrow(dist_df))
for (i in 1:nrow(dist_df)) {
  list_hldr[[i]] <- lapply(
    rand_walk_metrics$thinDraws,
    FUN = ks.test,
    y = dist_df[i, 'dist'],
    shape = dist_df[i, 'shape'],
    rate = dist_df[i, 'rate']
  )
}


# saving the power test
pdf(file = "../../images_slice_sampler_comp/curve3_rand_walk_power.pdf")
extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
dev.off()

# printing out metrics table
saveRDS(rand_walk_metrics,paste0("../../data/curve",curve_num,"_rand_walk_metrics"))

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "../../images_slice_sampler_comp/curve3_rand_walk.pdf")
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


rand_walk_min_max <- results(rand_walk_metrics, method = "Random Walk")

