rm(list=ls())
# setwd("~/cucumber/sam_comparison")
library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)
source("sam_comparison_functions.R")


reps <- 1

fexp <- function(x) exp(lf(x))
samples <- rep(500000,reps)
auto.cor.lim <- 0.05


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


# max min values
results <- function(df, method){
  if(method == 'Stepping Out'){
    # max_param <- df[df$SampPSec == max(df$SampPSec),c('x','w')][1,]
    # min_param <- df[df$SampPSec == min(df$SampPSec),c('x','w')][1,]
    #
    # max <- df %>%
    #   filter(x == max_param$x & w == max_param$w) %>%
    #   group_by(x,w) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'Stepping Out')
    #
    # min <- df %>%
    #   filter(x == min_param$x & w == min_param$w) %>%
    #   group_by(x,w) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'Stepping Out')
    #
    # max_df <- data.frame('Tuning Parameters' = c('x,w'),
    #                      'Values' = paste(max$x,',',max$w),
    #                      'Samples' = max$Samples,
    #                      'Accepted' = max$Accepted)
    #
    # min_df <- data.frame('Tuning Parameters' = c('x,w'),
    #                      'Values' = paste(min$x,',',min$w),
    #                      'Samples' = min$Samples,
    #                      'Accepted' = min$Accepted)
    #
    # return(rbind(max_df, min_df) %>%
    #   separate('Values',into = c('x','w'), sep = ',') %>%
    #   mutate(x = paste('$x =',x,'$'),
    #          w = paste('$w =',w,'$')) %>%
    #   unite('Settings',x:w, sep = ',') %>%
    #   mutate('Method' = 'Stepping Out') %>%
    #   relocate('Method', .before = Settings) %>%
    #   select(-Tuning.Parameters))

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
    # max_param <- df[df$SampPSec == max(df$SampPSec),c('x','s','rate')][1,]
    # min_param <- df[df$SampPSec == min(df$SampPSec),c('x','s','rate')][1,]
    #
    # max <- df %>%
    #   filter(x == max_param$x & s == max_param$s & rate == max_param$rate) %>%
    #   group_by(x,s,rate) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'Latent')
    #
    #
    # min <- df %>%
    #   filter(x == min_param$x & s == min_param$s & rate == min_param$rate) %>%
    #   group_by(x,s,rate) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'Latent')
    #
    # max_df <- data.frame('Tuning Parameters' = c('x,s,rate'),
    #            'Values' = paste(max$x,',',max$s,',',max$rate),
    #            'Samples' = max$Samples,
    #            'Accepted' = max$Accepted)
    #
    # min_df <- data.frame('Tuning Parameters' = c('x,s,rate'),
    #            'Values' = paste(min$x,',',min$s,',',min$rate),
    #            'Samples' = min$Samples,
    #            'Accepted' = min$Accepted)
    #
    # return(rbind(max_df, min_df) %>%
    #          separate('Values',into = c('x','s','rate'), sep = ',') %>%
    #          mutate(x = paste('$x =',x,'$'),
    #                 s = paste('$s =',s,'$'),
    #                 rate = paste('$\\text{rate} =',rate,'$')) %>%
    #          unite('Settings',x:rate, sep = ',') %>%
    #          mutate('Method' = 'Latent') %>%
    #          relocate('Method', .before = Settings) %>%
    #          select(-Tuning.Parameters))

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
    # max_param <- df[df$SampPSec == max(df$SampPSec),c('x','mu','sigma','df')][1,]
    # min_param <- df[df$SampPSec == min(df$SampPSec),c('x','mu','sigma','df')][1,]
    #
    # max <- df %>%
    #   filter(x == max_param$x & mu == max_param$mu & sigma == max_param$sigma & df == max_param$df) %>%
    #   group_by(x,mu,sigma,df) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'GESS')
    #
    #
    # min <- df %>%
    #   filter(x == min_param$x & mu == min_param$mu & sigma == min_param$sigma & df == min_param$df) %>%
    #   group_by(x,mu,sigma,df) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'GESS')
    #
    # max_df <- data.frame('Tuning Parameters' = c('x,mu,sigma,df'),
    #                      'Values' = paste(max$x,',',max$mu,',',max$sigma,',',max$df),
    #                      'Samples' = max$Samples,
    #                      'Accepted' = max$Accepted)
    #
    # min_df <- data.frame('Tuning Parameters' = c('x,mu,sigma,df'),
    #                      'Values' = paste(min$x,',',min$mu,',',min$sigma,',',min$df),
    #                      'Samples' = min$Samples,
    #                      'Accepted' = min$Accepted)
    #
    # return(rbind(max_df, min_df) %>%
    #          separate('Values',into = c('x','mu','sigma','df'), sep = ',') %>%
    #          mutate(x = paste('$x =',x,'$'),
    #                 mu = paste('$\\mu =',mu,'$'),
    #                 sigma = paste('$\\sigma =',sigma,'$'),
    #                 df = paste('$\\text{df} =',df,'$')) %>%
    #          unite('Settings',x:df, sep = ',') %>%
    #          mutate('Method' = 'GESS') %>%
    #          relocate('Method', .before = Settings) %>%
    #          select(-Tuning.Parameters))

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
    # max_param <- df[df$SampPSec == max(df$SampPSec),c('x','inv_cdf')][1,]
    # min_param <- df[df$SampPSec == min(df$SampPSec),c('x','inv_cdf')][1,]
    #
    # max <- df %>%
    #   filter(x == max_param$x & inv_cdf == max_param$inv_cdf) %>%
    #   group_by(x,inv_cdf) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'Transform')
    #
    #
    # min <- df %>%
    #   filter(x == min_param$x & inv_cdf == min_param$inv_cdf) %>%
    #   group_by(x,inv_cdf) %>%
    #   summarise(Samples = mean(SampPSec),
    #             Accepted = mean(ksTest > 0.05)) %>%
    #   mutate(Method = 'Transform')
    #
    # max_df <- data.frame('Tuning Parameters' = c('x,inv_cdf'),
    #                      'Values' = paste(max$x,',',max$inv_cdf),
    #                      'Samples' = max$Samples,
    #                      'Accepted' = max$Accepted)
    #
    # min_df <- data.frame('Tuning Parameters' = c('x,inv_cdf'),
    #                      'Values' = paste(min$x,',',min$inv_cdf),
    #                      'Samples' = min$Samples,
    #                      'Accepted' = min$Accepted)
    #
    # return(rbind(max_df, min_df) %>%
    #          separate('Values',into = c('x','inv_cdf'), sep = ' ,') %>%
    #          mutate(x = paste('$x =',x,'$'),
    #                 inv_cdf = paste('$\\text{inv cdf} =$',inv_cdf)) %>%
    #          unite('Settings',x:inv_cdf, sep = ',') %>%
    #          mutate('Method' = 'Transform') %>%
    #          relocate('Method', .before = Settings) %>%
    #          select(-c(Tuning.Parameters)))

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



# dist_comp <- function(inv_cdf,support){
#   if(is.list(support)) {
#   temp <- sub('q','d',inv_cdf)
#   temp <- sub('u,','support[[1]],', temp)
#   px <- eval(parse(text = temp))
#   px
#   } else {
#   temp <- sub('q','d',inv_cdf)
#   temp <- sub('u,','support,', temp)
#   px <- eval(parse(text = temp))
#   px
#   }
# }

dist_comp <- function(inv_cdf, find_grid){
  support <- eval(parse(text = find_grid))
  temp <- sub('q','d',inv_cdf)
  temp <- sub('u,','support,', temp)
  px <- eval(parse(text = temp))
  px
}


find_support <- function(inv_cdf){
  temp <- sub('u,','1e-05/2,',inv_cdf)
  lwr <- eval(parse(text = temp))
  temp <- sub('u,','1e-05/2,',inv_cdf)
  temp <- sub(')',',lower.tail = FALSE)',temp)
  upr <- eval(parse(text = temp))
  support <- seq(from = lwr, to = upr, length.out = 1000)
  support
}
