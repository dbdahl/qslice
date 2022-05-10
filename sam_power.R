
# 3
############ Curve 3 ###############
## dgamma(x, shape=2.5, log=TRUE) ##

# curve 3 in unimodal skewed right
lf <- function(x) {
  dgamma(x, shape=2.5, log=TRUE)
}
# plot of the log of the density function below

curve(fexp(x), xlim = c(0,60), ylim = c(0,.31))


samples <- 10000000
w <- 2
x <- 2#c(0.5, 2, 5)


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
         thin = min(which(acf(draws)$acf <0.01)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws)) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)

shape_values <- seq(from = 2.5 - 2.4, to = 2.5 + 100, by = 0.1)
pval_stepping_out <- vector(length = length(shape_values))

for(i in 1:length(pval_stepping_out)) {
  pval_stepping_out[i] <- ks.test(stepping_out_metrics$thinDraws[[1]],
                      pgamma(x, shape = shape_values[i]))$p.value
}

pval_stepping_out <- unlist(pval_stepping_out)

data.frame(pvalue = pval_stepping_out, shape = shape_values) %>% 
  ggplot(aes(x = shape, y = pvalue)) +
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
         thin = min(which(acf(draws, lag.max = 1000)$acf <0.01)),
         thinDraws = list(LaplacesDemon::Thin(draws, thin)),
         samplesThin = length(thinDraws)) %>%
  select(-metrics) %>%
  mutate(SampPSec = ESS/time) %>%
  relocate(samples, .after = time)


pval_latent <- vector(length = length(shape_values))

for(i in 1:length(pval_latent)) {
  pval_latent[i] <- ks.test(latent_metrics$thinDraws[[1]],
                      pgamma(x, shape = shape_values[i]))$p.value
}

pval_latent <- unlist(pval_latent)

data.frame(pvalue = pval_latent, shape = shape_values) %>% 
  ggplot(aes(x = shape, y = pvalue)) +
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


pval_gess <- vector(length = length(shape_values))

for(i in 1:length(pval_gess)) {
  pval_gess[i] <- ks.test(gess_metrics$thinDraws[[1]],
                      pgamma(x, shape = shape_values[i]))$p.value
}

pval_gess <- unlist(pval_gess)

data.frame(pvalue = pval_gess, shape = shape_values) %>% 
  ggplot(aes(x = shape, y = pvalue)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = 'red')


data.frame(shape_values, pval_stepping_out, pval_latent, pval_gess) %>% 
  reshape2::melt(id.vars = 'shape_values', variable.name = 'sampler', value.name = 'pval') %>% 
  ggplot(aes(x = shape_values, y = pval, col = sampler)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = 'red') + 
  geom_vline(xintercept = 2.5, color = 'blue')
