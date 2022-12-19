rm(list=ls())
# setwd("~/cucumber/sam_comparison")
library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)
source("../sam_comparison_functions.R")


reps <- 10

fexp <- function(f, x) exp(f(x))
samples <- rep(5000,reps)
auto.cor.lim <- 0.10
sampleSize <- 10000

# function that extracts pvalues when calculating power
extract_pvals <- function(list_hldr, dist_df) {
  counter <- 1
  pvalue_list <- list(length = nrow(dist_df))
  temp_df <- data.frame(matrix(nrow = nrow(dist_df) * reps, ncol = ncol(dist_df)))
  for(j in 1:nrow(dist_df)) {
    pvalue <- vector(length = reps)
    for(i in 1:reps) {
      pvalue[[i]] <- list_hldr[[j]][[i]]$p.value
      temp_df[counter,] <- dist_df[j,]
      counter <- counter + 1
    }
    pvalue_list[[j]] <- pvalue
  }
  cbind(temp_df, pVal = unlist(pvalue_list)) %>%
    unite('Distribution', X1:X3) %>%
    ggplot(aes(x = pVal, color = Distribution)) +
    geom_histogram() +
    facet_wrap(~Distribution) +
    theme(
      legend.position = 'none'
    )
}


# returns the results of any slice sampling method in a nice table
results <- function(df, method){
  if(method == 'Stepping Out'){
    param_settings <- df %>%
      group_by(x,w) %>%
      summarise(Samples = mean(SampPSec),
                Accepted = mean(ksTest > 0.05)) %>%
      mutate(Method = 'Stepping Out')

    param_settings_df <- data.frame('Tuning Parameters' = c('x,w'),
                         'Values' = paste(param_settings$x,',',param_settings$w),
                         'Samples' = param_settings$Samples,
                         'Accepted' = param_settings$Accepted)

    return(param_settings_df %>%
      separate('Values',into = c('x','w'), sep = ',') %>%
      mutate(x = paste('$x =',x,'$'),
             w = paste('$w =',w,'$')) %>%
      unite('Settings',x:w, sep = ',') %>%
      mutate('Method' = 'Stepping Out') %>%
      relocate('Method', .before = Settings) %>%
      select(-Tuning.Parameters))
  }
  if(method == 'Latent'){
    param_settings <- df %>%
      group_by(x,s,rate) %>%
      summarise(Samples = mean(SampPSec),
                Accepted = mean(ksTest > 0.05)) %>%
      mutate(Method = 'Latent')
    
    param_settings_df <- data.frame('Tuning Parameters' = c('x,s,rate'),
               'Values' = paste(param_settings$x,',',param_settings$s,',',param_settings$rate),
               'Samples' = param_settings$Samples,
               'Accepted' = param_settings$Accepted)
    
    return(param_settings_df %>%
             separate('Values',into = c('x','s','rate'), sep = ',') %>%
             mutate(x = paste('$x =',x,'$'),
                    s = paste('$s =',s,'$'),
                    rate = paste('$\\text{rate} =',rate,'$')) %>%
             unite('Settings',x:rate, sep = ',') %>%
             mutate('Method' = 'Latent') %>%
             relocate('Method', .before = Settings) %>%
             select(-Tuning.Parameters))
  }
  if(method == 'GESS'){
    param_settings <- df %>%
      group_by(x,mu,sigma,df) %>%
      summarise(Samples = mean(SampPSec),
                Accepted = mean(ksTest > 0.05)) %>%
      mutate(Method = 'GESS')

    param_settings_df <- data.frame('Tuning Parameters' = c('x,mu,sigma,df'),
                         'Values' = paste(param_settings$x,',',param_settings$mu,',',param_settings$sigma,',',param_settings$df),
                         'Samples' = param_settings$Samples,
                         'Accepted' = param_settings$Accepted)

    return(param_settings_df%>%
             separate('Values',into = c('x','mu','sigma','df'), sep = ',') %>%
             mutate(x = paste('$x =',x,'$'),
                    mu = paste('$\\mu =',mu,'$'),
                    sigma = paste('$\\sigma =',sigma,'$'),
                    df = paste('$\\text{df} =',df,'$')) %>%
             unite('Settings',x:df, sep = ',') %>%
             mutate('Method' = 'GESS') %>%
             relocate('Method', .before = Settings) %>%
             select(-Tuning.Parameters))
  }
  if(method == 'Transform'){
    param_settings <- df %>%
      group_by(x,inv_cdf) %>%
      summarise(Samples = mean(SampPSec),
                Accepted = mean(ksTest > 0.05)) %>%
      mutate(Method = 'Transform')

    param_settings_df <- data.frame('Tuning Parameters' = c('x,inv_cdf'),
                         'Values' = paste(param_settings$x,',',param_settings$inv_cdf),
                         'Samples' = param_settings$Samples,
                         'Accepted' = param_settings$Accepted)

    return(param_settings_df %>%
             separate('Values',into = c('x','inv_cdf'), sep = ' ,') %>%
             mutate(x = paste('$x =',x,'$'),
                    inv_cdf = paste('$\\text{inv cdf} =$',inv_cdf)) %>%
             unite('Settings',x:inv_cdf, sep = ',') %>%
             mutate('Method' = 'Transform') %>%
             relocate('Method', .before = Settings) %>%
             select(-c(Tuning.Parameters)))
  }
  if(method == 'Random Walk') {
    param_settings <- df %>%
      group_by(x,c) %>%
      summarise(Samples = mean(SampPSec),
                Accepted = mean(ksTest > 0.05)) %>%
      mutate(Method = 'Random Walk')

    param_settings_df <- data.frame('Tuning Parameters' = c('x,c'),
                                    'Values' = paste(param_settings$x,',',param_settings$c),
                                    'Samples' = param_settings$Samples,
                                    'Accepted' = param_settings$Accepted)

    return(param_settings_df %>%
             separate('Values',into = c('x','c'), sep = ' ,') %>%
             mutate(x = paste('$x =',x,'$'),
                    c = paste('$c =$',c)) %>%
             unite('Settings',x:c, sep = ',') %>%
             mutate('Method' = 'Random Walk') %>%
             relocate('Method', .before = Settings) %>%
             select(-c(Tuning.Parameters)))
  }
}

# this function returns the results of the transform slice sampler in a nice table
results_transform <- function(df) {
  param_settings <- df %>%
    group_by(x,inv_cdf) %>%
    summarise(KLD = mean(KLD),
              JSD = mean(JSD),
              Samples = mean(SampPSec),
              Accepted = mean(ksTest > 0.05)) %>%
    mutate(Method = 'Transform')

  param_settings_df <- data.frame('Tuning Parameters' = c('x,inv_cdf'),
                                  'Values' = paste(param_settings$x,',',param_settings$inv_cdf),
                                  'KLD' = param_settings$KLD,
                                  'JSD' = param_settings$JSD,
                                  'Samples' = param_settings$Samples,
                                  'Accepted' = param_settings$Accepted)

  return(param_settings_df %>%
           separate('Values',into = c('x','inv_cdf'), sep = ' ,') %>%
           mutate(x = paste('$x =',x,'$'),
                  inv_cdf = paste('$\\text{inv cdf} =$',inv_cdf)) %>%
           unite('Settings',x:inv_cdf, sep = ',') %>%
           mutate('Method' = 'Transform') %>%
           relocate('Method', .before = Settings) %>%
           select(-c(Tuning.Parameters)))
}

# this function finds the second derivative of function f at point x
second_derivative <- function( x, h = 1e-5, f ) {
  num <- f(x + h) - 2*f(x) + f(x - h)
  denom <- h^2
  num/denom
}

# this function provides a lapace approximation for a curve lf
laplace_approx <- function(lf, h = 1e-5, init){
  fit <- optim(init, lf, control=list(fnscale=-1), method = 'BFGS')
  sd <- sqrt(-solve(second_derivative(fit$par, h = h, f = lf))) * 3
  return(
      c(log_pdf = function(x) {dnorm(x, mean = fit$par, sd = sd, log = TRUE)},
      inv_cdf = function(u, lower.tail = TRUE) {qnorm(u, mean = fit$par, sd = sd, lower.tail = lower.tail)})
  )
}

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
