rm(list=ls())

library(magrittr)
library(tidyverse)
library(cucumber)
# source("functions_log.R")
source("sam_comparison_functions.R")

auto.cor.lim <- 0.05

options(xtable.format.args = list(big.mark = ","), # separates large numbers using a ,
        # xtable.size = "\\tiny", # makes the font tiny to fit into the overleaf presentation
        xtable.append = FALSE, # replaces the table
        xtable.table.placement = 'h', # table placment h is for "here"
        xtable.caption.placement = 'bottom', # places the caption at the bottom
        xtable.include.rownames = FALSE, # makes it so there are now row names
        # xtable.hline.after = seq(from = -1, to = nrow(tab), by = 1), # adds a line after every row
        xtable.sanitize.text.function = function(x){x}, # makes it so latex shows in the column names
        adjustcolor.alpha.f = 0.02 #adjustcolor('black', alpha.f = 0.95) # changes the opacity of the lines comparing target to samples
        )


fexp <- function(x) exp(lf(x))
samples <- 1000000

####
## stepping out metrics to input ##
w <- c(2)#c(0.01, 1, 2, 4, 10)

##
## latent slice sampling metric to input ##
s <- c(3)#c(0.01, 1, 2, 10)
rate <- c(2)#c(0.5, 1, 1.5, 2, 2.5, 3)

####
## gess slice sampling metrics to input ##
mu <- c(3)#c(1,2,3,4.5,6,7)
sigma <- c(3)#c(2,3,4,5,6,8)
df <- c(3)#c(1,4,16,16^2,16^4)

##
## elliptical slice sampling metrics to input ##
# there are none besides starting point


# 1
################### Curve 1 ########################
## 0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2) ##


# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2))
}

cdf <- function(x) {
  probs <- vector(length = length(x))
  for(i in 1:length(x)){
    probs[i] <- integrate(fexp, lower = -Inf, upper = x[i])$value
  }
  return(probs)
}

# plot of the log of the density function below
pdf(file = "images_slice_sampler_comp/curve1.pdf")
curve(fexp(x), xlim = c(-4,15), ylim = c(0,.17))
dev.off()

xlim_range <- c(-4,15)
ylim_range <- c(0,0.17 + 0.10)

x <- c(5,0)#c(0.00001, 2, 5, 7, 10, 15)
##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out metrics table
tab <- xtable::xtable(stepping_out_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                       'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve1_metrics_tbl_stepping_out.tex',
                     # ,
                     size = "\\tiny")

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve1_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve1_w_stepping_out.tex')


# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve1_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black'))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws))%>% rename('n' = 'samples',
                                                                                'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve1_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve1_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve1_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve1_rate_latent.tex')


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve1_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##


# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                    x_0 = x,
                                    mu_value = mu,
                                    sigma_value = sigma,
                                    df_value = df,
                                    lf_func = lf,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out metrics table
tab <- xtable::xtable(gess_metrics %>% select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                               "$\\sigma$" = 'sigma',
                                                                               'n' = 'samples',
                                                                               'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve1_metrics_tbl_gess.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve1_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve1_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve1_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve1_df_gess.tex')




# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve1_gess.pdf")
curve(fexp(x),col = 'red',xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red',xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

# 2
################## Curve2 ###########################
## 0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=20,sd=1) ##

# curve 2 is bimodal
lf <- function(x) {
  log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=20,sd=1))
}

cdf <- function(x) {
  probs <- vector(length = length(x))
  for(i in 1:length(x)){
    probs[i] <- integrate(fexp, lower = -Inf, upper = x[i])$value
  }
  return(probs)
}
# plot of the log of the density function below
pdf(file = "images_slice_sampler_comp/curve2.pdf")
curve(fexp(x), xlim = c(-3,25), ylim = c(0,.32))
dev.off()

x <- c(0, 20)

xlim_range <- c(-3,25)
ylim_range <- c(0,0.32 + 0.10)

##
## Stepping out
##

# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(stepping_out_metrics %>% select(-c(draws, thinDraws))%>% rename('n' = 'samples',
                                                                                      'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve2_metrics_tbl_stepping_out.tex',
                     size = "\\tiny")


# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve2_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve2_w_stepping_out.tex')

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve2_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                 'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve2_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve2_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve2_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve2_rate_latent.tex')


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve2_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 1000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(gess_metrics %>% select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                               "$\\sigma$" = 'sigma',
                                                                               'n' = 'samples',
                                                                               'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve2_metrics_tbl_gess.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve2_start_gess.tex',
      )

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve2_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve2_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve2_df_gess.tex')


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve2_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

# 3
############ Curve 3 ###############
## dgamma(x, shape=2.5, log=TRUE) ##

# curve 3 in unimodal skewed right
lf <- function(x) {
  dgamma(x, shape=2.5, rate = 1, log=TRUE)
}
# plot of the log of the density function below
pdf(file = "images_slice_sampler_comp/curve3.pdf")
curve(fexp(x), xlim = c(0,10), ylim = c(0,.31))
dev.off()

x <- 2#c(0.5, 2, 5)

xlim_range <- c(0,10)
ylim_range <- c(0,0.31 + 0.10)

##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pgamma, shape = 2.5, rate = 1)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)


# printing out the metrics table
tab <- xtable::xtable(stepping_out_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                       'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve3_metrics_tbl_stepping_out.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve3_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve3_w_stepping_out.tex')

pdf(file = "images_slice_sampler_comp/curve3_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.99))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##


# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 1000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pgamma, shape = 2.5, rate = 1)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                 'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve3_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve3_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve3_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve3_rate_latent.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve3_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.99))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##


# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 1000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pgamma, shape = 2.5, rate = 1)$p.value) %>%
  dplyr::select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(gess_metrics %>% dplyr::select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                                      "$\\sigma$" = 'sigma',
                                                                                      'n' = 'samples',
                                                                                      'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve3_metrics_tbl_gess.tex',
                     size = "\\tiny")


# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve3_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve3_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve3_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve3_df_gess.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve3_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

# 4
######################### Curve 4 ############################
## ifelse( x < 0, -Inf, dt(x, df=3.0, log=TRUE) + log(2.0)) ##

# curve 4 is right skewed with support for x strictly positive
lf <- function(x) {
  ifelse( x < 0, -Inf, dt(x, df=3.0, log=TRUE) + log(2.0))
}

cdf <- function(x) {
  probs <- vector(length = length(x))
  for(i in 1:length(x)){
    probs[i] <- integrate(fexp, lower = 0, upper = x[i])$value
  }
  return(probs)
}

# plot of the log of the density function below
pdf(file = "images_slice_sampler_comp/curve4.pdf")
curve(fexp(x), xlim = c(0,10), ylim = c(0,.8))
dev.off()

x <- 1#c(0.5,1,2,4)

xlim_range <- c(0,10)
ylim_range <- c(0,0.8 + 0.10)

##
## Stepping out
##

# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  dplyr::select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(stepping_out_metrics %>% dplyr::select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                       'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve4_metrics_tbl_stepping_out.tex',
                     size = "\\tiny")



# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve4_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve4_w_stepping_out.tex')


pdf(file = 'images_slice_sampler_comp/curve4_stepping_out.pdf')
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                 'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve4_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve4_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve4_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve4_rate_latent.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve4_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  lwr_support = 0,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, cdf)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(gess_metrics %>% select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                               "$\\sigma$" = 'sigma',
                                                                               'n' = 'samples',
                                                                               'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve4_metrics_tbl_gess.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve4_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve4_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve4_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve4_df_gess.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve4_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()


beepr::beep()

# 5
########## Curve 5 ##########
## dt(x, df=3.0, log=TRUE) ##

# t distribution with a steep peak
lf <- function(x) {
  dt(x, df=3.0, log=TRUE)
}
# plot of the log of the density function below
pdf(file = "images_slice_sampler_comp/curve5.pdf")
curve(fexp(x), xlim = c(-10,10), ylim = c(0,.42))
dev.off()

x <- 0#c(0, 2, 5)

xlim_range <- c(-10,10)
ylim_range <- c(0,0.42 + 0.10)


##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pt, df = 3.0)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(stepping_out_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                       'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve5_metrics_tbl_stepping_out.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve5_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve5_w_stepping_out.tex')

# evalTbl_stepping_out <- cbind(start_points_stepping_out, w_values_stepping_out) %>% round(.,1)

pdf(file = "images_slice_sampler_comp/curve5_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##


# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pt, df = 3.0)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                 'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve5_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve5_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve5_s_latent.tex')
# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve5_rate_latent.tex')


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve5_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##


# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pt, df = 3.0)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(gess_metrics %>% select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                               "$\\sigma$" = 'sigma',
                                                                               'n' = 'samples',
                                                                               'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve5_metrics_tbl_gess.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve5_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve5_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve5_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve5_df_gess.tex')


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve5_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()


# 6
################### Curve 6 ####################
## dbeta(x, shape1=0.2, shape2=0.8, log=TRUE) ##

# u shaped distribution
lf <- function(x) {
  dbeta(x, shape1=0.2, shape2=0.8, log=TRUE)
}
# plot of the log of the density function below
pdf(file = "images_slice_sampler_comp/curve6.pdf")
curve(fexp(x), xlim = c(0,1), ylim = c(0,10))
dev.off()

x <- 0.2#c(0.1, 0.4, 0.6, 0.9)

xlim_range <- c(0,1)
ylim_range <- c(0, 5 + 0.10)


##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 1000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pbeta, shape1 = 0.2, shape2 = 0.8)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(stepping_out_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                       'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve6_metrics_tbl_stepping_out.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve6_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve6_w_stepping_out.tex')

pdf(file = "images_slice_sampler_comp/curve6_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pbeta, shape1 = 0.2, shape2 = 0.8)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                 'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve6_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve6_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve6_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve6_rate_latent.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve6_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  upr_support = 1,
                                  lwr_support = 0,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws, lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pbeta, shape1 = 0.2, shape2 = 0.8)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(gess_metrics %>% select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                               "$\\sigma$" = 'sigma',
                                                                               'n' = 'samples',
                                                                               'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve6_metrics_tbl_gess.tex',
                     size = "\\tiny")


# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve6_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve6_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve6_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve6_df_gess.tex')


pdf(file = "images_slice_sampler_comp/curve6_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

beepr::beep(sound = 8)


# 7
############## Curve 7 ######################
## dnorm(x, mean = 20, sd = 5, log = TRUE) ##

# normal prior normal likelihood function
lf <- function(x) {
  dnorm(x, mean = 20, sd = 5, log = TRUE)
}
# plot of the density !! can use for elliptical slice sampler !!
pdf(file = "images_slice_sampler_comp/curve7.pdf")
curve(fexp(x), xlim = c(0,40))
dev.off()

x <- 20#c(10, 20, 30)

xlim_range <- c(0,40)
ylim_range <- c(0, 0.08 + 0.010)

##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pnorm, mean = 20, sd = 5)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(stepping_out_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                       'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve7_metrics_tbl_stepping_out.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve7_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve7_w_stepping_out.tex')

pdf(file = "images_slice_sampler_comp/curve7_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws,lag.max = 5000)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pnorm, mean = 20, sd = 5)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(latent_metrics %>% select(-c(draws, thinDraws)) %>% rename('n' = 'samples',
                                                                                 'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve7_metrics_tbl_latent.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve7_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve7_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve7_rate_latent.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve7_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pnorm, mean = 20, sd = 5)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(gess_metrics %>% select(-c(draws, thinDraws)) %>% rename("$\\mu$" = 'mu',
                                                                               "$\\sigma$" = 'sigma',
                                                                               'n' = 'samples',
                                                                               'nThin' = 'samplesThin'),
                      digits = c(0,0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve7_metrics_tbl_gess.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve7_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve7_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve7_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve7_df_gess.tex')


pdf(file = "images_slice_sampler_comp/curve7_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

###
### elliptical
###

lf <- function(x) {
  emdbook::dmvnorm(x,
                   mu = t(matrix(c(10,
                                 5), nrow = 2, ncol = 1, byrow = TRUE)),
                   Sigma = matrix(c(25,0,
                                    0,25), nrow = 2, ncol = 2, byrow = TRUE),
                   log = TRUE)
}


cdf <- function(x) {
  probs <- vector(length = length(x))
  for(i in 1:length(x)){
    probs[i] <- integrate(fexp, lower = -Inf, upper = x[i])$value
  }
  return(probs)
}


slice_sampler_elliptical <- function(x, target, mu=2, sigma=5, log=TRUE) {
  if ( ! isTRUE(log) ) stop("'log=FALSE' is not implemented.")
  nEvaluations <- 0
  f <- function(x) { nEvaluations <<- nEvaluations + 1; target(x) }
  # Step 1
  y <- log(runif(1)) + f(x)
  nu <- MASS::mvrnorm(n = 1, mu = mu, Sigma = sigma)#rnorm(1,mu,sigma)
  theta <- runif(1,0,2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  repeat {
    x1 <- (x - mu)*cos(theta) + (nu - mu)*sin(theta) + mu
    if ( y < f(x1) ) return(list(x=x1, nEvaluations=nEvaluations))
    if (theta < 0) theta_min <- theta
    else theta_max <- theta
    theta <- runif(1, theta_min, theta_max)
  }
}

# these are the values for the normal prior distribution. Only used for the elliptical slice sampler
mu_prior <- matrix(c(5,
                     0), nrow = 2, ncol = 1, byrow = TRUE)
sigma_prior <- matrix(c(4,0,
                        0,4), nrow = 2, ncol = 2, byrow = TRUE)
x_start_vector <- matrix(c(2,
                           0), nrow = 2, ncol = 1, byrow = TRUE)

samp1 <- slice_sampler_elliptical(x = t(x_start_vector),
                         target = lf,
                         mu = t(mu_prior),
                         sigma = sigma_prior,
                         log = TRUE)

samp2 <- slice_sampler_elliptical(x = samp1$x,
                                  target = lf,
                                  mu = t(mu_prior),
                                  sigma = sigma_prior,
                                  log = TRUE)


elliptical_samples <- elliptical_time_eval(samples = samples,
                                           lf_func = lf,
                                           x_0 = x_start_vector,
                                           mu_value = t(mu_prior),
                                           sigma_value = sigma_prior,
                                           log_value = TRUE) %>%
  rowwise() %>%
  mutate(thinDraws = list(LaplacesDemon::Thin(Draws, 5)))


# creating a data frame with all possible combinations
trials_elliptical <- expand.grid(samples, x_start_vector, mu_prior, sigma_prior) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4')



# creating a data frame for the metrics
elliptical_metrics <- trials_elliptical %>%
  dplyr::rowwise() %>%
  mutate(metrics = elliptical_time_eval(samples = samples,
                                        lf_func = lf,
                                        x_0 = x,
                                        mu_value = mu,
                                        sigma_value = sigma,
                                        log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         ksTest = ks.test(thinDraws, pnorm, mean = 20, sd = 5)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# printing out the metrics table
tab <- xtable::xtable(elliptical_metrics %>% select(-c(draws, thinDraws)),
                      digits = c(0,0,0,0,0,0,2,0,0,0,2,0))
xtable::print.xtable(tab, file = 'images_slice_sampler_comp/curve7_metrics_tbl_elliptical.tex',
                     size = "\\tiny")

# evaluation of each starting point
start_points_elliptical <- elliptical_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_elliptical)
print(tab, file = 'images_slice_sampler_comp/curve7_start_elliptical.tex')


pdf(file = "images_slice_sampler_comp/curve7_elliptical.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(elliptical_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

beepr::beep()


# 8
############################# Curve 8 ################################
## log(dnorm(x, mean = 4, sd = 2) * dgamma(x, shape = 2, rate = 4)) ##

# normal prior with gamma likelihood function
lf <- function(x) {
  log(dnorm(x, mean = 4, sd = 2) * dgamma(x, shape = 2, rate = 4))
}
# plot of the density curve !! can use for elliptical slice sampler !!
pdf(file = "images_slice_sampler_comp/curve8.pdf")
curve(fexp(x), xlim = c(-5,10))
dev.off()

xlim_range <- c(-5,10)
ylim_range <- c(0, 0.05 + 0.010)

x <- c(0.00001, 1, 2)

# these are the values for the normal prior distribution. Only used for the elliptical slice sampler
mu_prior <- 4
sigma_prior <- 2

##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve8_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve8_w_stepping_out.tex')

pdf(file = "images_slice_sampler_comp/curve8_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve8_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve8_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve8_rate_latent.tex')



# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve8_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve8_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\mu$" = 'mu')

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve8_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec)) %>%
  rename("$\\sigma$" = 'sigma')

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve8_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve8_df_gess.tex')


pdf(file = "images_slice_sampler_comp/curve8_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

###
### elliptical
###

# creating a data frame with all possible combinations
trials_elliptical <- expand.grid(samples, x, mu_prior, sigma_prior) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4')

# creating a data frame for the metrics
elliptical_metrics <- trials_elliptical %>%
  dplyr::rowwise() %>%
  mutate(metrics = elliptical_time_eval(samples = samples,
                                        lf_func = lf,
                                        x_0 = x,
                                        mu_value = mu,
                                        sigma_value = sigma,
                                        log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_elliptical <- elliptical_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve8_elliptical.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(elliptical_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

# 9
################################## Curve 9 ###################################
## log(dnorm(x, mean = 4, sd = 2) * (dbeta(x, shape1 = .10, shape2 = .05))) ##

# normal prior with beta likelihood function
lf <- function(x) {
  log(dnorm(x, mean = 4, sd = 2) * (dbeta(x, shape1 = .10, shape2 = .05)))
}
# plot of the density curve !! can use for elliptical slice sampler !!
pdf(file = 'images_slice_sampler_comp/curve9.pdf')
curve(fexp(x), xlim = c(0,1), ylim = c(0,.25))
dev.off()

xlim_range <- c(0,1)
ylim_range <- c(0, 0.15 + 0.10)

x <- c(0.00001, 0.2, 0.6, 0.8, 0.9999)

# these are the values for the normal prior distribution. Only used for the elliptical slice sampler
mu_prior <- 4
sigma_prior <- 2

##
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3')

# evaluating the sampler at each different iteration
stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = lf,
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_stepping_out <- stepping_out_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve9_start_stepping_out.tex')

# evaluation at each value of w
w_values_stepping_out <- stepping_out_metrics %>%
  group_by(w) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(w_values_stepping_out)
print(tab, file = 'images_slice_sampler_comp/curve9_w_stepping_out.tex')

pdf(file = "images_slice_sampler_comp/curve9_stepping_out.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(stepping_out_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## latent
##

# creating a data frame with all possible combinations
trials_latent <- expand.grid(samples, x, s, rate) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                's' = 'Var3',
                'rate' = 'Var4')

# creating a data frame for the metrics
latent_metrics <- trials_latent %>%
  dplyr::rowwise() %>%
  mutate(metrics = latent_time_eval(samples = samples,
                                    x_0 = x,
                                    s_0 = s,
                                    lf_func = lf,
                                    rate_value = rate,
                                    log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_latent <- latent_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_latent)
print(tab, file = 'images_slice_sampler_comp/curve9_start_latent.tex')

# evaluation at each value of s
s_values_latent <- latent_metrics %>%
  group_by(s) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(s_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve9_s_latent.tex')

# evaluation at each rate value
rate_values_latent <- latent_metrics %>%
  group_by(rate) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(rate_values_latent)
print(tab, file = 'images_slice_sampler_comp/curve9_rate_latent.tex')


pdf(file = "images_slice_sampler_comp/curve9_latent.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(latent_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

##
## gess
##

# creating a data frame with all possible combinations
trials_gess <- expand.grid(samples, x, mu, sigma, df) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4',
                'df' = 'Var5')

# creating a data frame for the metrics
gess_metrics <- trials_gess %>%
  dplyr::rowwise() %>%
  mutate(metrics = gess_time_eval(samples = samples,
                                  x_0 = x,
                                  mu_value = mu,
                                  sigma_value = sigma,
                                  df_value = df,
                                  lf_func = lf,
                                  log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_gess <- gess_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(start_points_gess)
print(tab, file = 'images_slice_sampler_comp/curve9_start_gess.tex')

# evaluation at each value of mu
mu_values_gess <- gess_metrics %>%
  group_by(mu) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(mu_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve9_mu_gess.tex')

# evaluation at each sigma value
sigma_values_gess <- gess_metrics %>%
  group_by(sigma) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(sigma_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve9_sigma_gess.tex')

# evaluation at each df value
df_values_gess <- gess_metrics %>%
  group_by(df) %>%
  summarise(avgSampPSec = mean(SampPSec))

tab <- xtable::xtable(df_values_gess)
print(tab, file = 'images_slice_sampler_comp/curve9_df_gess.tex')


pdf(file = "images_slice_sampler_comp/curve9_gess.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(gess_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()

###
### elliptical
###

# creating a data frame with all possible combinations
trials_elliptical <- expand.grid(samples, x, mu_prior, sigma_prior) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'mu' = 'Var3',
                'sigma' = 'Var4')

# creating a data frame for the metrics
elliptical_metrics <- trials_elliptical %>%
  dplyr::rowwise() %>%
  mutate(metrics = elliptical_time_eval(samples = samples,
                                        lf_func = lf,
                                        x_0 = x,
                                        mu_value = mu,
                                        sigma_value = sigma,
                                        log_value = TRUE))  %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf < auto.cor.lim)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin))) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

# evaluation of each starting point
start_points_elliptical <- elliptical_metrics %>%
  group_by(x) %>%
  summarise(avgSampPSec = mean(SampPSec))


# evalTbl_latent <- cbind(start_points_latent, s_values_latent, rate_values_latent) %>% round(.,1)
pdf(file = "images_slice_sampler_comp/curve9_elliptical.pdf")
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2)
lapply(elliptical_metrics$thinDraws, function(x) {
  lines(density(x), col = adjustcolor('black', alpha.f = 0.95))
})
curve(fexp(x),col = 'red', xlim = c(xlim_range[1],xlim_range[2]), ylim = c(ylim_range[1],ylim_range[2]), lwd = 2, add = TRUE)
dev.off()


####
############ Creating a Loop for each function ################
####


fexp <- function(x) exp(lf(x))

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2))
}
# plot of the log of the density function below
curve(lf(x), xlim = c(-15,15))

# setting values to evaluate the stepping out procedure
lf_functions <- c('log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2))',
                  'log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=20,sd=1))',
                  'dgamma(x, shape=2.5, log=TRUE)',
                  'ifelse( x < 0, -Inf, dt(x, df=1.0, log=TRUE) + log(2.0))',
                  'dt(x, df=1.0, log=TRUE)',
                  'dbeta(x, shape1=0.2, shape2=0.8, log=TRUE)',
                  'log(dnorm(x, mean = 1, sd = 1) * dnorm(x, mean = 5, sd = 9))',
                  'log(dnorm(x, mean = 4, sd = 2) * dgamma(x, shape = 2, rate = 4))',
                  'log(dnorm(x, mean = 4, sd = 2) * (dbeta(x, shape1 = .10, shape2 = .05)))'
)
samples <- 500
x <- c(0.01, 0.2, 0.5, 0.8, 0.99)
w <- c(0.01, 1, 2, 10, 100)

# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w, lf_functions) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3',
                'lf' = 'Var4')

stepping_out_metrics <- trials_stepping_out %>%
  rowwise() %>%
  mutate(metrics = stepping_out_time_eval(samples = samples,
                                          x_0 = x,
                                          lf_func = function(x){eval(parse(text = lf))},
                                          w_value = w,
                                          max_value = Inf,
                                          log_value = TRUE)) %>%
  mutate(nEval = metrics$nEval,
         ESS = metrics$EffSamp,
         time = metrics$Time) %>%
  select(-metrics) %>%
  relocate(samples, .after = time)


stepping_out_time_eval(samples = 500,
                       x_0 = 0.4,
                       lf_func = function(x){eval(parse(text=lf_functions))} ,
                       w_value = 4,
                       max_value = Inf,
                       log_value = TRUE)

