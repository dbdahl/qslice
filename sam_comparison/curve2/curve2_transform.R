# setwd('~/cucumber/sam_comparison/curve2')
source('curve2_setup.R')

##
#### Transform ####
##

px <- data.frame(px = matrix(nrow = length(log_pdf), ncol = 1)) 
px$log_pdf <- log_pdf
px$grid <- find_grid
px <- px %>%
  rowwise() %>% 
  mutate(px = list(exp(log_pdf(grid)))) %$% px

temp_df <- data.frame(log_pdf = matrix(nrow = length(log_pdf), ncol = 1))
temp_df$log_pdf <- log_pdf
temp_df$inv_cdf <- inv_cdf
temp_df$px <- px

# creating a data frame with all possible combinations
trials_transform <- expand.grid(samples, x) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2')

trials_transform <-
  sapply(trials_transform, rep.int, times = length(log_pdf)) %>% data.frame()

transform_parameters <-
  sapply(temp_df, rep.int, times = nrow(trials_transform) / length(log_pdf))  %>%
  data.frame() %>%
  arrange(inv_cdf)

trials_transform <- cbind(trials_transform, transform_parameters)

# evaluating the sampler at each different iteration
transform_metrics <- trials_transform %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    metrics = transform_time_eval(
      samples = samples,
      x_0 = x,
      lf_func = lf,
      pseudo_pdf_log = log_pdf,
      pseudo_cdf_inv = inv_cdf,
      log_value = TRUE
    )
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
    ksTest = ks.test(truncThinDraws, cdf)$p.value,
    KLD.JSD = list(LaplacesDemon::KLD(px = px, py = py)[c(4,6)]),
    KLD = KLD.JSD$sum.KLD.px.py,
    JSD = KLD.JSD$mean.sum.KLD
  ) %>%
  dplyr::select(-c(metrics,KLD.JSD)) %>%
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
#     transform_metrics$thinDraws,
#     FUN = ks.test,
#     y = cdf,
#     mean2 = dist_df[i, 'mean2'],
#     sd2 = dist_df[i, 'sd2']
#   )
# }
# 
# 
# # saving the power test
# pdf(file = "../../images_slice_sampler_comp/curve2_transform_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(transform_metrics,paste0("../../data/curve",curve_num,"_transform_metrics"))

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "../../images_slice_sampler_comp/curve2_transform.pdf")
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2
)
lapply(transform_metrics$thinDraws, function(x) {
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


transform_min_max <-
  results(transform_metrics, method = "Transform")

transform_evaluation <- results_transform(transform_metrics)

tab <- transform_evaluation %>%
  xtable::xtable(caption = "This table shows the method of sampling, settings for the tuning parameters, average number of effective samples taken per second, the perecentage of samples that matched the target disribution according to a kolmogorov-Smirnov test, and the Kullback-Leibler divergence.",
                 digits = c(0,0,0,2,2,0,2))

xtable::print.xtable(tab, file = '../../images_slice_sampler_comp/curve2_transform_evaluation.tex')


rm(
  # list_hldr,
  # dist_df,
  trials_transform,
  transform_metrics
)