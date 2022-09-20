# setwd("~/cucumber/sam_comparison")
source('setup.R')

lf <- function(x) {
  dnorm(x, mean = 4, sd = 4)
}

cdf <- function(x) {
  pnorm(x, 4, 4)
}

log_pdf <- c('dnorm(x, mean = 5, sd = 5, log = TRUE)',
             'dnorm(x, mean = 4, sd = 10, log = TRUE)')
inv_cdf <- c('qnorm(u, mean = 5, sd = 5)',
             'qnorm(u, mean = 4, sd = 10)')
samples <- rep(1000, 10)
x <- c(4, 2)
w <- c(1, 4)
s <- c(2, 6)
rate <- c(2)
mu <- c(0, 4)
sigma <- c(1, 4)
df <- c(3)
method <- c('Stepping Out', 'Latent', 'GESS', 'Transform')

max(method %in% 'Stepping Out')


sam <- slice_sampler_comparison(
  lf = lf,
  cdf = cdf,
  curve_num = 10,
  samples = samples,
  x = x,
  w = w,
  s = s,
  rate = rate,
  mu = mu,
  sigma = sigma,
  df = df,
  inv_cdf = inv_cdf,
  log_pdf = log_pdf,
  method = method
)

slice_sampler_comparison <- function(
           lf, # needs to be the log the target dist function
           cdf, # the cdf of the true value
           curve_num, # what you want to call the curve
           samples, # number of samples
           x , # starting points
           w = NULL, # w values
           s = NULL, # s values
           rate = NULL, # rate values
           mu = NULL, # mu values
           sigma = NULL, # sigma values
           df = NULL, # df values
           inv_cdf = NULL, # inv cdf
           log_pdf = NULL, # lot of pseudo prior
           method){
    # plot of the log of the density function below
    pdf(file = paste0("../images_slice_sampler_comp/curve", curve_num, ".pdf"))
    curve(fexp(x))
    dev.off()

    ##
    #### Stepping out ####
    ##

    if (max(method %in% 'Stepping Out') == 1) {
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

      # printing out metrics table
      tab <-
        xtable::xtable(
          stepping_out_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                           'nThin' = 'samplesThin'),
          digits = c(0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0)
        )
      xtable::print.xtable(
        tab,
        file = paste0(
          '../images_slice_sampler_comp/curve',
          curve_num,
          '_metrics_tbl_stepping_out.tex'
        ),
        size = "\\tiny"
      )


      # evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
      pdf(file = paste0(
        "../images_slice_sampler_comp/curve",
        curve_num,
        "_stepping_out.pdf"
      ))
      curve(fexp(x), col = 'red', lwd = 2)
      lapply(stepping_out_metrics$thinDraws, function(x) {
        lines(density(x), col = adjustcolor('black'))
      })
      curve(fexp(x),
            col = 'red',
            lwd = 2,
            add = TRUE)
      dev.off()

      stepping_min_max <-
        results(stepping_out_metrics, method = "Stepping Out")

      rm(trials_stepping_out,
         stepping_out_metrics)

      print('Finished Stepping Out')

    }


    ##
    #### latent ####
    ##
    if (max(method %in% 'Latent') == 1) {
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


      # printing out metrics table
      tab <-
        xtable::xtable(
          latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                     'nThin' = 'samplesThin'),
          digits = c(0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0)
        )
      xtable::print.xtable(
        tab,
        file = paste0(
          '../images_slice_sampler_comp/curve',
          curve_num,
          '_metrics_tbl_latent.tex'
        ),
        size = "\\tiny"
      )


      pdf(file = paste0(
        "../images_slice_sampler_comp/curve",
        curve_num,
        "_latent.pdf"
      ))
      curve(fexp(x), col = 'red', lwd = 2)
      lapply(latent_metrics$thinDraws, function(x) {
        lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
      })
      curve(fexp(x),
            col = 'red',
            lwd = 2,
            add = TRUE)
      dev.off()


      latent_min_max <- results(latent_metrics, method = "Latent")

      rm(trials_latent,
         latent_metrics)

      print('Finished Latent')

    }

    ##
    #### gess ####
    ##

    if (max(method %in% 'GESS') == 1) {
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


      # printing out metrics table
      tab <-
        xtable::xtable(
          gess_metrics %>% select(-c(draws, thinDraws)) %>% rename(
            "$\\mu$" = 'mu',
            "$\\sigma$" = 'sigma',
            'n' = 'samples',
            'nThin' = 'samplesThin'
          ),
          digits = c(0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0)
        )
      xtable::print.xtable(
        tab,
        file = paste0(
          '../images_slice_sampler_comp/curve',
          curve_num,
          '_metrics_tbl_gess.tex'
        ),
        size = "\\tiny"
      )


      # evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
      pdf(file = paste0(
        "../images_slice_sampler_comp/curve",
        curve_num,
        "_gess.pdf"
      ))
      curve(fexp(x), col = 'red', lwd = 2)
      lapply(gess_metrics$thinDraws, function(x) {
        lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
      })
      curve(fexp(x),
            col = 'red',
            lwd = 2,
            add = TRUE)
      dev.off()

      gess_min_max <- results(gess_metrics, method = "GESS")


      rm(trials_gess,
         gess_metrics)

      print('Finished GESS')

    }
    ##
    #### Transform ####
    ##

    if (max(method %in% 'Transform') == 1) {


      temp_df <- data.frame(log_pdf,
                            inv_cdf)

      # creating a data frame with all possible combinations
      trials_transform <- expand.grid(samples, x) %>%
        dplyr::rename('samples' = 'Var1',
                      'x' = 'Var2')

      trials_transform <-
        sapply(trials_transform, rep.int, times = length(log_pdf)) %>% data.frame()

      transform_parameters <-
        sapply(temp_df,
               rep.int,
               times = nrow(trials_transform) / length(log_pdf))

      trials_transform <- cbind(trials_transform, transform_parameters)

      # evaluating the sampler at each different iteration
      transform_metrics <- trials_transform %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          metrics = transform_time_eval(
            samples = samples,
            x_0 = x,
            lf_func = lf,
            pseudo_pdf_log = function(x) {
              eval(parse(text = log_pdf))
            },
            pseudo_cdf_inv = function(u) {
              eval(parse(text = inv_cdf))
            },
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
          support = find_support(inv_cdf),
          KLD.JSD = list(LaplacesDemon::KLD(px = dist_comp(inv_cdf = inv_cdf, support = support), py = py)[c(4,6)]),
          KLD = KLD.JSD$sum.KLD.px.py,
          JSD = KLD.JSD$mean.sum.KLD
        ) %>%
        dplyr::select(-c(metrics, KLD.JSD)) %>%
        dplyr::mutate(SampPSec = ESS / time) %>%
        dplyr::relocate(samples, .after = time)


      # printing out metrics table
      tab <-
        xtable::xtable(
          transform_metrics %>% select(-c(draws, thinDraws))  %>% rename('n' = 'samples',
                                                                         'nThin' = 'samplesThin'),
          digits = c(0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0)
        )
      xtable::print.xtable(
        tab,
        file = paste0(
          '../images_slice_sampler_comp/curve',
          curve_num,
          '_metrics_tbl_transform.tex'
        ),
        size = "\\tiny"
      )


      # evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
      pdf(file = paste0(
        "../images_slice_sampler_comp/curve",
        curve_num,
        "_transform.pdf"
      ))
      curve(fexp(x), col = 'red', lwd = 2)
      lapply(transform_metrics$thinDraws, function(x) {
        lines(density(x), col = adjustcolor('black'))
      })
      curve(fexp(x),
            col = 'red',
            lwd = 2,
            add = TRUE)
      dev.off()

      transform_min_max <-
        results(transform_metrics, method = "Transform")

      rm(trials_transform,
         transform_metrics)

      print('Finished Transform')

    }

    # creating min max table
    min_max <- rbind(stepping_min_max,
                     latent_min_max,
                     gess_min_max,
                     transform_min_max)

    # printing out min max table
    tab <- min_max %>%
      xtable::xtable(caption = 'This table shows the method of sampling, settings for the tuning parameters, average number of effective samples taken per second and the perecentage of samples that matched the target disribution according to a kolmogorov-Smirnov test. This table shows that for some of the sampling methods, the effectivness of the method is largely dependon the choice of tunning parameters.',
                     digits = c(0, 0, 0, 0, 2))

    xtable::print.xtable(tab,
                         file = paste0(
                           '../images_slice_sampler_comp/curve',
                           curve_num,
                           '_min_max.tex'
                         ))

    min_max

  }
