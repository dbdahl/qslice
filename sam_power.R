library(tidyverse)
source("sam_comparison_functions.R")

fexp <- function(x) exp(lf(x))

# 3
############ Curve 3 ###############
## dgamma(x, shape=2.5, log=TRUE) ##

# curve 3 in unimodal skewed right
lf <- function(x) {
  #dgamma(x, shape=2.5, log=TRUE)
  dt(x, df=20.0, log=TRUE)
}

# cdf <- function(x) {
#   probs <- vector(length = length(x))
#   for(i in 1:length(x)){
#     probs[i] <- integrate(fexp, lower = -Inf, upper = x[i])$value
#   }
#   return(probs)
# }

# plot of the density function below
curve(fexp(x), xlim = c(-3,3))

x <- c(1)

parameter_values <- c(1,5,10,15,20,25,30,35,40)

samples <- rep(100000,100)#c(1000000,10000,1000,100)

# setting values to evaluate the stepping out procedure
lf_functions <- c(paste0("dt(x, df=", parameter_values, ", log = TRUE)"),
                  "dnorm(x, mean=0,sd=1, log = TRUE)",
                  "dnorm(x, mean=2,sd=1, log = TRUE)",
                  "dgamma(x, shape = 6, rate = 3, log = TRUE)")
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
## Stepping out
##


# creating a data frame with all possible combinations
trials_stepping_out <- expand.grid(samples, x, w, lf_functions) %>%
  dplyr::rename('samples' = 'Var1',
                'x' = 'Var2',
                'w' = 'Var3',
                'lf' = 'Var4') %>%
  mutate(lf = as.character(lf))

# evaluating the sampler at each different iteration
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
         time = metrics$Time,
         draws = metrics$Draws,
         thin = min(which(acf(draws)$acf <0.01)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws),
         pVals = ks.test(x = thinDraws, y = pt, df = 20.0)$p.value) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

pdf(file = 'images_slice_sampler_comp/power_histogram.pdf')
stepping_out_metrics %>%
  ggplot(aes(x = pVals, color = as.factor(lf))) +
  geom_histogram() +
  facet_wrap(~as.factor(lf)) +
  geom_vline(xintercept = 0.05, color = 'red') +
  theme(
    legend.position = 'none'
  )
dev.off()



samp_num <- rep(100,100)

rand_func <- c(paste0("rt(n = samples, df=", parameter_values, ")"),
                  "rnorm(n = samples, mean=0,sd=1)",
                  "rnorm(n = samples, mean=2,sd=1)",
                  "rgamma(n = samples, shape = 6, rate = 3)")


rand_df <- expand.grid(samp_num, rand_func) %>%
  dplyr::rename('samples' = 'Var1',
                'lf' = 'Var2') %>%
  mutate(lf = as.character(lf)) %>%
  rowwise() %>%
  mutate(draws = list(eval(parse(text = gsub("samples", samples, lf))))) %>%
  mutate(pVals = ks.test(draws, pt, df = 20)$p.value)


rand_df %>%
  ggplot(aes(x = pVals, color = as.factor(lf))) +
  geom_histogram() +
  facet_wrap(~as.factor(lf)) +
  geom_vline(xintercept = 0.05, color = 'red') +
  theme(
    legend.position = 'none'
  )

#####non-finished#####

power <- list()

for(i in parameter_values){
  pvals <- sapply(stepping_out_metrics$thinDraws,
                function(data) {ks.test(x = data, y = pt(x, df = i))$p.value})
  power[[i]] <- pvals
}

pval_stepping_out <- vector(length = length(parameter_values))

for(i in 1:length(pval_stepping_out)) {
  pval_stepping_out[i] <- ks.test(stepping_out_metrics$thinDraws[[1]],
                      pt(x, df = parameter_values[i]))$p.value
}

pval_stepping_out <- unlist(pval_stepping_out)

data.frame(pvalue = pval_stepping_out, parameter = parameter_values) %>%
  ggplot(aes(x = parameter, y = pvalue)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = 'red')



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
         thin = min(which(acf(draws, lag.max = 2000)$acf <0.01)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws)) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)


pval_latent <- vector(length = length(parameter_values))

for(i in 1:length(pval_latent)) {
  pval_latent[i] <- ks.test(latent_metrics$thinDraws[[1]],
                      pt(x, df = parameter_values[i]))$p.value
}

pval_latent <- unlist(pval_latent)

data.frame(pvalue = pval_latent, parameter = parameter_values) %>%
  ggplot(aes(x = parameter, y = pvalue)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = 'red')


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
         thin = min(which(acf(draws, lag.max = 1000)$acf <0.01)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws)) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)


pval_gess <- vector(length = length(parameter_values))

for(i in 1:length(pval_gess)) {
  pval_gess[i] <- ks.test(gess_metrics$thinDraws[[1]],
                      pt, df = parameter_values[i])$p.value
}

pval_gess <- unlist(pval_gess)

data.frame(pvalue = pval_gess, parameter = parameter_values) %>%
  ggplot(aes(x = parameter, y = pvalue)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = 'red')



###
###
###

pdf(file = "images_slice_sampler_comp/power_curve_gamma.pdf")
data.frame(parameter_values, pval_stepping_out, pval_latent, pval_gess) %>%
  reshape2::melt(id.vars = 'parameter_values', variable.name = 'sampler', value.name = 'pval') %>%
  ggplot(aes(x = parameter_values, y = pval, col = sampler)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = 'red') +
  geom_vline(xintercept = 3.0, color = 'blue') +
  geom_hline(yintercept = ks.test(stepping_out_metrics$thinDraws[[1]],
                                  pnorm, mean = 0, sd = 1)$p.value, color = 'gold', size = 1.2) +
  geom_hline(yintercept = ks.test(stepping_out_metrics$thinDraws[[1]],
                                  pbeta, shape1 = 0.2, shape2 = 0.8)$p.value, color = 'darkgreen', size = 1.2) +
  geom_hline(yintercept = ks.test(stepping_out_metrics$thinDraws[[1]],
                                  pt, df = 3)$p.value, color = 'purple', size = 1.2) +
  scale_x_continuous(limits = c(0,15))
dev.off()




  ###
  ###
  ###

# pval_stepping_out_samples <- data.frame(matrix(nrow = 520, ncol = 3)) %>%
#   rename('samples' = 'X1',
#          'pval' = 'X2',
#          'parameter' = 'X3')
# row <- 1
#
# for(i in 1:length(parameter_values)) {
#   for(j in 1:4) {
#     pval_stepping_out_samples[row, 1] <- length(stepping_out_metrics$thinDraws[[j]])
#     pval_stepping_out_samples[row, 2] <- ks.test(stepping_out_metrics$thinDraws[[j]],
#                                     pt(x, df = parameter_values[i]))$p.value
#     pval_stepping_out_samples[row, 3] <- parameter_values[i]
#     row <- row + 1
#   }
# }


pval_stepping_out_samples <- data.frame(matrix(nrow = 520, ncol = 3)) %>%
  rename('samples' = 'X1',
         'pval' = 'X2',
         'parameter' = 'X3')
sampSize <- c(300000,100000,10000,1000)
row <- 1

for(i in 1:length(parameter_values)) {
  for(j in 1:4) {
    pval_stepping_out_samples[row, 1] <- sampSize[j]
    pval_stepping_out_samples[row, 2] <- ks.test(sample(stepping_out_metrics$thinDraws[[1]], sampSize[j]),
                                                 pt, df = parameter_values[i])$p.value
    pval_stepping_out_samples[row, 3] <- parameter_values[i]
    row <- row + 1
  }
}

pdf(file = 'images_slice_sampler_comp/power_curve_samples_t.pdf')
pval_stepping_out_samples %>%
  ggplot(aes(x = parameter, y = pval, color = as.factor(samples))) +
  geom_point() +
  geom_vline(xintercept = 3.0, color = 'blue')+
  labs(
    color = 'samples'
  ) +
  geom_hline(yintercept = 0.05, color = 'red')
dev.off()


pval_stepping_out <- vector(length = length(parameter_values))
for(j in 1:length(sampSize)) {
  for(i in 1:length(pval_stepping_out)) {
    pval_stepping_out[i] <- ks.test(sample(stepping_out_metrics$thinDraws[[1]], sampSize[j]),
                                    pt, df = parameter_values[i])$p.value
  }
}

pval_stepping_out <- unlist(pval_stepping_out)



pdf(file = 'images_slice_sampler_comp/gamma_distributions.pdf')
plot(density(stepping_out_metrics$thinDraws[[1]]), col = 'red', xlim = c(0,7), ylim = c(0,0.60))
curve(dgamma(x, shape = 3), add = TRUE, col = 'blue')
curve(dgamma(x, shape = 2), add = TRUE, col = 'green')
curve(dgamma(x, shape = 1), add = TRUE, col = 'brown')
curve(dgamma(x, shape = 0.1), add = TRUE, col = 'pink')
dev.off()
