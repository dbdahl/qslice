#### loading in library ####
library(tidyverse)
library(stringr)
library(patchwork)
library(rlist)
library(ggthemes)

### theme ###
theme_set(
  theme_classic()+
    theme(
      panel.grid = element_blank(),
      legend.position = 'none'
    )
)


# my_gg_colors <- RColorBrewer::brewer.pal(6, 'Dark2') 
# options(ggplot2.discrete.colour = my_gg_colors,
#         ggplot2.discrete.fill = my_gg_colors)

# function for eda #
fexp <- function(f, x) exp(f(x))

extract_curve_metrics <- function(path) {
  temp <- readRDS(path)
  temp <- temp %>% 
    mutate(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
                          method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk')) %>% 
    select(-c(draws,thinDraws))
  list(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
       method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk'),
       data = temp)
}

boxplot_slice_sampler <- function(list_Of_df, curve_num) {
  
  stepping_out <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'stepping_out')
  stepping_out <- stepping_out[[1]]$data
  
  # stepping out
  stepping_out_x <- stepping_out %>% 
    ggplot(aes(y = as.factor(x), x = ESS)) +
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  stepping_out_w <- stepping_out %>% 
    ggplot(aes(y = as.factor(w), x = ESS)) + 
    geom_boxplot() + 
    labs(
      x = 'Effective Sample Size',
      y = 'w'
    )+
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  stepping_out_param <- stepping_out %>% 
    mutate(parameters = paste('x =',x,'w =',w)) %>% 
    ggplot(aes(y = as.factor(parameters), x = ESS)) + 
    geom_boxplot() + 
    labs(
      x = 'Effective Sample Size',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
    
  # latent
  
  latent <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'latent')
  latent <- latent[[1]]$data

  latent_x <- latent %>% 
    ggplot(aes(y = as.factor(x), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  latent_s <- latent %>% 
    ggplot(aes(y = as.factor(s), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 's'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  latent_rate <- latent %>% 
    ggplot(aes(y = as.factor(rate), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'rate'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  latent_param <- latent %>% 
    mutate(parameters = paste('x = ',x,'s =',s,'rate =',rate)) %>% 
    ggplot(aes(x = ESS, y = as.factor(parameters))) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  # gess
  gess <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'gess')
  gess <- gess[[1]]$data

  gess_x <- gess %>% 
    ggplot(aes(y = as.factor(x), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  gess_mu <- gess %>% 
    ggplot(aes(y = as.factor(mu), x = ESS)) +
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'mu'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  gess_sigma <- gess %>% 
    ggplot(aes(y = as.factor(sigma), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'sigma'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  gess_df <- gess %>% 
    ggplot(aes(y = as.factor(df), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'df'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  gess_param <- gess %>% 
    mutate(parameters = paste('x =',x,'mu =',mu,'sigma =',sigma,'df =',df)) %>% 
    ggplot(aes(x = ESS, y = as.factor(parameters))) +
    geom_boxplot() + 
    labs(
      x = 'Effective Sample Size',
      y = 'parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  # transform
  transform <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'transform')
  transform <- transform[[1]]$data

  
  transform$inv_cdf_chr <- as.character(transform$inv_cdf) %>% str_extract('\\bq.*')
  
  transform_x <- transform %>% 
    ggplot(aes(y = as.factor(x), x = ESS)) + 
    geom_boxplot()  +
    labs(
      x = 'Effective Sample Size',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  transform_inv_cdf <- transform %>% 
    ggplot(aes(y = as.factor(inv_cdf_chr), x = ESS, fill = as.factor(KLD))) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'Inverse CDF'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  transform_param <- transform %>% 
    mutate(parameters = paste('x =',x,'inv_cdf =',inv_cdf_chr)) %>% 
    ggplot(aes(y = as.factor(parameters), x = ESS)) +
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  # random walk
  rand_walk <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'rand_walk')
  rand_walk <- rand_walk[[1]]$data
  
  rand_walk_x <- rand_walk %>% 
    ggplot(aes(y = as.factor(x), x = ESS)) +
    geom_boxplot()  +
    labs(
      x = 'Effective Sample Size',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  rand_walk_c <- rand_walk %>% 
    ggplot(aes(y = as.factor(c), x = ESS)) +
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'c'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  rand_walk_param <- rand_walk %>% 
    mutate(parameters = paste('x =',x,'c =',c)) %>% 
    ggplot(aes(y = as.factor(parameters), x = ESS)) + 
    geom_boxplot() +
    labs(
      x = 'Effective Sample Size',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  
  
  # total comp
  sampler_comp <- data.frame(
    method = c(rep('stepping_out', nrow(stepping_out)),
               rep('latent',nrow(latent)),
               rep('gess', nrow(gess)),
               rep('transform', nrow(transform)),
               rep('rand_walk', nrow(rand_walk))),
    ESS = c(stepping_out$ESS, latent$ESS, gess$ESS, transform$ESS, rand_walk$ESS)
  ) %>% 
    ggplot(aes(y = as.factor(method), x = ESS)) +
    geom_boxplot()  +
    labs(
      x = 'Effective Sample Size',
      y = 'Method'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K'))
  

  return(list(
    stepping_out_x = stepping_out_x, stepping_out_w = stepping_out_w, stepping_out_param = stepping_out_param,
    stepping_out_plot = ((stepping_out_x + stepping_out_w)/stepping_out_param),
    latent_x = latent_x, latent_s = latent_s, latent_rate = latent_rate, latent_param = latent_param,
    latent_plot = ((latent_x + latent_s)/(latent_rate)/(latent_param)),
    gess_x = gess_x, gess_mu = gess_mu, gess_sigma = gess_sigma, gess_df = gess_df, gess_param = gess_param,
    gess_plot = ((gess_x + gess_mu)/(gess_sigma + gess_df)/(gess_param)),
    transform_x = transform_x, transform_inv_cdf = transform_inv_cdf, transform_param = transform_param,
    transform_plot = ((transform_x)/ (transform_inv_cdf)/(transform_param)),
    rand_walk_x = rand_walk_x, rand_walk_c = rand_walk_c, rand_walk_param = rand_walk_param, sampler_comp = sampler_comp,
    rand_walk_plot = ((rand_walk_x + rand_walk_c)/(rand_walk_param))
  ))
}


extract_metrics <- function(path) {
  temp <- readRDS(path)
  temp <- temp %>% mutate(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
                          method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk')) %>% 
    select(method, curve, ESS)
}



extract_draws <- function(path) {
  temp <- readRDS(path)
  temp <- temp %>% 
    mutate(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
           method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk')) %>% 
    select(curve, method, draws, thinDraws)
  list(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
       method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk'),
       data = temp)
}


curve_slice_sampler <- function(list_Of_df, curve) {
  
  curves <- c(
    function(x) {
      log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 6, sd = 2))
    },
    function(x) {
      log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 20, sd = 1))
    },
    function(x) {
      dgamma(x,
             shape = 2.5,
             rate = 1,
             log = TRUE)
    },
    function(x) {
      ifelse(x < 0, -Inf, dt(x, df = 3.0, log = TRUE) + log(2.0))
    },
    function(x) {
      dt(x, df = 5.0, log = TRUE)
    },
    function(x) {
      dbeta(x,
            shape1 = 0.2,
            shape2 = 0.8,
            log = TRUE)
    },
    function(x) {
      dnorm(x, mean = 20, sd = 5, log = TRUE)
    },
    function(x) {
      if_else(x < 0, log(0),
              if_else(x <= 5, log((x^3 + x^2 + 100*sin(x))/(-100*cos(5) + (3575)/(12))),
                      log(0)))
    }
  )
  
  range <- list( 
    list(xrange = c(-4, 15), yrange = c(0, 0.17 + 0.10)),
    list(xrange = c(-3, 25), yrange = c(0, 0.32 + 0.10)),
    list(xrange = c(0, 10), yrange = c(0, 0.31 + 0.10)),
    list(xrange = c(0, 10), yrange = c(0, 0.8 + 0.10)),
    list(xrange = c(-10, 10), yrange = c(0, 0.42 + 0.10)),
    list(xrange = c(0, 1), yrange = c(0, 5 + 0.10)),
    list(xrange = c(0, 40), yrange = c(0, 0.08 + 0.010)),
    list(xrange = c(0,5), yrange = c(0, 0.4))   
  )
  
  lf <- curves[[curve]]
  
  par(mfrow = c(2,3))
  # stepping out
  stepping_out <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'stepping_out') 
  stepping_out <- stepping_out[[curve]]$data
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2
  )
  lapply(stepping_out$thinDraws, function(x) {
    lines(density(x), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2,
    add = TRUE
  )
  
  # latent 
  latent <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'latent')
  latent <- latent[[curve]]$data
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2
  )
  lapply(latent$thinDraws, function(x) {
    lines(density(x), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2,
    add = TRUE
  )
  
  # gess
  gess <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'gess')
  gess <- gess[[curve]]$data
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2
  )
  lapply(gess$thinDraws, function(x) {
    lines(density(x), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2,
    add = TRUE
  )
  
  # transform
  transform <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'transform')
  transform <- transform[[curve]]$data
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2
  )
  lapply(transform$thinDraws, function(x) {
    lines(density(x), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2,
    add = TRUE
  )
  
  # rand_walk
  rand_walk <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'rand_walk')
  rand_walk <- rand_walk[[curve]]$data
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2
  )
  lapply(rand_walk$thinDraws, function(x) {
    lines(density(x), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    lwd = 2,
    add = TRUE
  )
}

