# setwd("~/cucumber/sam_comparison")
source('setup.R')
curve_num <- 8

# 8
################### Curve 8 #############################
## (x^3 + x^2 + 100*sin(x))(-100*cos(5) + (3575)/(12)) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  if_else(x < 0, log(0),
          if_else(x <= 5, log((x^3 + x^2 + 100*sin(x))/(-100*cos(5) + (3575)/(12))),
                  log(0)))
}

cdf <- function(x) {
  probs <- vector(length = length(x))
  for (i in 1:length(x)) {
    probs[i] <- if_else(x <0, 0,
                        if_else(x <= 5, ((12)/(-1200*cos(5) + 3575))*((x^4)/(4) + (x^3)/(3) + 100*(-cos(x) + 1)), 1))
  }
  return(probs)
}

# plot of the log of the density function below
pdf(file = "../images_slice_sampler_comp/curve8.pdf")
curve(fexp(f = lf, x = x), xlim = c(0,5), ylim = c(0, 0.4))
dev.off()

grid <- seq(from = 0,
            to = 5,
            length.out = 1000)

py <- exp(lf(grid))

xlim_range <- c(0,3)
ylim_range <- c(0, 0.4)

#### Tuning Parameters ####
## starting point ##
x <- c(0.5, 1, 4.5)
## stepping out metrics to input ##
w <- c(2, 5, 10)

## latent slice sampling metric to input ##
s <- c(3, 5, 10)
rate <- c(2)

## gess slice sampling metrics to input ##
mu <- c(3, 5, 10)
sigma <- c(3)
df <- c(3)

## random walk tuning parameters ##
c <- c(1)

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
  dplyr::rowwise() %>%
  dplyr::mutate(
    metrics = stepping_out_time_eval(
      samples = samples,
      x_0 = x,
      lf_func = lf,
      w_value = w,
      max_value = Inf,
      log_value = TRUE
    )
  ) %>%
  dplyr::mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = min(which(
      acf(draws, plot = FALSE, lag.max = 1000)$acf < auto.cor.lim
    )),
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, cdf)$p.value
  ) %>%
  dplyr::select(-metrics) %>%
  dplyr::mutate(SampPSec = ESS / time) %>%
  dplyr::relocate(samples, .after = time)

# # the functions to compare against
# dist_df <- data.frame(
#   dist = rep("cdf", 6),
#   mean2 = c(6, 6, 6, 7, 5, 6),
#   sd2 = c(4, 1, 3, 2, 2, 2)
# )
# 
# list_hldr <- list(length = nrow(dist_df))
# for (i in 1:nrow(dist_df)) {
#   list_hldr[[i]] <- lapply(
#     stepping_out_metrics$thinDraws,
#     FUN = ks.test,
#     y = cdf,
#     mean2 = dist_df[i, 'mean2'],
#     sd2 = dist_df[i, 'sd2']
#   )
# }
# 
# 
# # saving the power test
# pdf(file = "../images_slice_sampler_comp/curve8_stepping_out_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(stepping_out_metrics,paste0("../data/curve",curve_num,"_stepping_out_metrics"))

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "../images_slice_sampler_comp/curve8_stepping_out.pdf")
curve(
  fexp(f = lf, x = x),
  col = 'red',
  xlim = c(xlim_range[1], xlim_range[2]),
  ylim = c(ylim_range[1], ylim_range[2]),
  lwd = 2
)
lapply(stepping_out_metrics$thinDraws, function(x) {
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

stepping_min_max <-
  results(stepping_out_metrics, method = "Stepping Out")

rm(
  list_hldr,
  dist_df,
  trials_stepping_out,
  stepping_out_metrics
)

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
    ksTest = ks.test(thinDraws, cdf)$p.value
  ) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# # the functions to compare against
# dist_df <- data.frame(
#   dist = rep("cdf", 6),
#   mean2 = c(6, 6, 6, 7, 5, 6),
#   sd2 = c(4, 1, 3, 2, 2, 2)
# )
# 
# list_hldr <- list(length = nrow(dist_df))
# for (i in 1:nrow(dist_df)) {
#   list_hldr[[i]] <- lapply(
#     latent_metrics$thinDraws,
#     FUN = ks.test,
#     y = cdf,
#     mean2 = dist_df[i, 'mean2'],
#     sd2 = dist_df[i, 'sd2']
#   )
# }
# 
# # saving the power test
# pdf(file = "../images_slice_sampler_comp/curve8_latent_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(latent_metrics,paste0("../data/curve",curve_num,"_latent_metrics"))

# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "../images_slice_sampler_comp/curve8_latent.pdf")
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
    ksTest = ks.test(thinDraws, cdf)$p.value
  ) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS / time) %>%
  relocate(samples, .after = time)

# # the functions to compare against
# dist_df <- data.frame(
#   dist = rep("cdf", 6),
#   mean2 = c(6, 6, 6, 7, 5, 6),
#   sd2 = c(4, 1, 3, 2, 2, 2)
# )
# 
# list_hldr <- list(length = nrow(dist_df))
# for (i in 1:nrow(dist_df)) {
#   list_hldr[[i]] <- lapply(
#     gess_metrics$thinDraws,
#     FUN = ks.test,
#     y = cdf,
#     mean2 = dist_df[i, 'mean2'],
#     sd2 = dist_df[i, 'sd2']
#   )
# }
# 
# # saving the power test
# pdf(file = "../images_slice_sampler_comp/curve8_gess_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(gess_metrics,paste0("../data/curve",curve_num,"_gess_metrics"))

# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "../images_slice_sampler_comp/curve8_gess.pdf")
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


##
#### Transform ####
##

log_pdf <- c(function(x) dnorm(x, mean = 5, sd = 5, log = TRUE),
             function(x) dnorm(x, mean = 4, sd = 10, log = TRUE))

inv_cdf <- c(function(u) qnorm(u, mean = 5, sd = 5),
             function(u) qnorm(u, mean = 4, sd = 10))

find_grid <- list(seq(from = qnorm((1-.99999)/2, mean = 5, sd = 5),
                      to = qnorm((1-.99999)/2, mean = 5, sd = 5, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 4, sd = 10),
                      to = qnorm((1-.99999)/2, mean = 4, sd = 10, lower.tail = FALSE), length.out = 1000))

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
  sapply(temp_df, rep.int, times = nrow(trials_transform) / length(log_pdf)) %>%
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
    thin = min(which(
      acf(draws, plot = FALSE, lag.max = 1000)$acf < auto.cor.lim
    )),
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, cdf)$p.value,
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
#   mean2 = c(6, 6, 6, 7, 5, 6),
#   sd2 = c(4, 1, 3, 2, 2, 2)
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
# pdf(file = "../images_slice_sampler_comp/curve8_transform_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(transform_metrics,paste0("../data/curve",curve_num,"_transform_metrics"))

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "../images_slice_sampler_comp/curve8_transform.pdf")
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

xtable::print.xtable(tab, file = '../images_slice_sampler_comp/curve8_transform_evaluation.tex')

rm(
  tab,
  list_hldr,
  dist_df,
  trials_transform,
  transform_metrics
)

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
    thin = min(which(
      acf(draws, plot = FALSE, lag.max = 5000)$acf < auto.cor.lim
    )),
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, cdf)$p.value
  ) %>%
  dplyr::select(-c(metrics)) %>%
  dplyr::mutate(SampPSec = ESS / time) %>%
  dplyr::relocate(samples, .after = time)

# # the functions to compare against
# dist_df <- data.frame(
#   dist = rep("cdf", 6),
#   mean2 = c(6, 6, 6, 7, 5, 6),
#   sd2 = c(4, 1, 3, 2, 2, 2)
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
# pdf(file = "../images_slice_sampler_comp/curve8_rand_walk_power.pdf")
# extract_pvals(list_hldr = list_hldr, dist_df = dist_df)
# dev.off()

# printing out metrics table
saveRDS(rand_walk_metrics,paste0("../data/curve",curve_num,"_rand_walk_metrics"))

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "../images_slice_sampler_comp/curve8_rand_walk.pdf")
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


# creating min max table
min_max <- rbind(stepping_min_max,
                 latent_min_max,
                 gess_min_max,
                 transform_min_max,
                 rand_walk_min_max)

# printing out min max table
tab <- min_max %>%
  xtable::xtable(caption = 'This table shows the method of sampling, settings for the tuning parameters, average number of effective samples taken per second and the perecentage of samples that matched the target disribution according to a kolmogorov-Smirnov test. This table shows that for some of the sampling methods, the effectivness of the method is largely dependon the choice of tunning parameters.',
                 digits = c(0, 0, 0, 0, 2))

xtable::print.xtable(tab, file = '../images_slice_sampler_comp/curve8_min_max.tex')
