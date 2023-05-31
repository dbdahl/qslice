## simulation study of the Cauchy transform
## author: Sam Johnson

# setwd('~/cucumber/sam_comparison/curve1')
source('norm_setup.R')

##
#### Transform ####
##


temp_df <- data.frame(log_pdf = matrix(nrow = length(log_pdf), ncol = 1))
# temp_df$scales <- c(scales,cauchyFit$sc)
temp_df$log_pdf <- log_pdf
temp_df$inv_cdf <- inv_cdf
temp_df$t <- t
temp_df$px <- lapply(temp_df$log_pdf, FUN = \(func) exp(func(grid)))

# creating a data frame with all possible combinations
trials_transform <- expand.grid(samples, x) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2')

trials_transform <-
  sapply(trials_transform, rep.int, times = nrow(temp_df)) %>% data.frame()

transform_parameters <-
  sapply(temp_df, rep.int, times = nrow(trials_transform) / nrow(temp_df)) %>%
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
  mutate(
    nEval = metrics$nEval,
    ESS = metrics$EffSamp,
    time = metrics$Time,
    draws = metrics$Draws,
    thin = (length(draws)/ESS) * 10,
    thinDraws = list(LaplacesDemon::Thin(draws, thin)),
    samplesThin = length(thinDraws),
    ksTest = ks.test(thinDraws, pnorm, 0, 1)$p.value,
    KLD.JSD = list(LaplacesDemon::KLD(px = px, py = py)[c(4,6)]),
    KLD = KLD.JSD$sum.KLD.px.py,
    JSD = KLD.JSD$mean.sum.KLD
  ) %>%
  dplyr::select(-c(metrics, KLD.JSD, thinDraws)) %>%
  dplyr::mutate(SampPSec = ESS / time) %>%
  dplyr::relocate(samples, .after = time)


transform_metrics <- transform_metrics %>% 
  mutate(time = time,#ifelse(grepl('*Auto*', t), time+burnin_metrics$Time, time),
         nEval = nEval#ifelse(grepl('*Auto*', t), nEval+burnin_metrics$nEval, nEval),
         )

if(saveInd) saveRDS(transform_metrics, file = 'data/transform.rds')
