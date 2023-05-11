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


# curves and ranges
curves <- c(
  function(x) {
    dnorm(x, mean = 20, sd = 5, log = TRUE)
  },
  function(x) {
    dgamma(x,
           shape = 2.5,
           rate = 1,
           log = TRUE)
  },
  function(x) {
    log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 6, sd = 2))
  },
  function(x) {
    log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 20, sd = 1))
  },
  function(x) {
    dt(x, df = 5.0, log = TRUE)
  },
  function(x) {
    ifelse(x < 0, -Inf, dt(x, df = 3.0, log = TRUE) + log(2.0))
  },
  function(x) {
    dbeta(x,
          shape1 = 0.2,
          shape2 = 0.8,
          log = TRUE)
  },
  function(x) {
    if_else(x < 0, log(0),
            if_else(x <= 5, log((x^3 + x^2 + 100*sin(x))/(-100*cos(5) + (3575)/(12))),
                    log(0)))
  }
)

range <- list( 
  list(xrange = c(0, 40), yrange = c(0, 0.08)),
  list(xrange = c(0, 10), yrange = c(0, 0.31)),
  list(xrange = c(-4, 15), yrange = c(0, 0.17)),
  list(xrange = c(-3, 25), yrange = c(0, 0.32)),
  list(xrange = c(-10, 10), yrange = c(0, 0.42)),
  list(xrange = c(0, 10), yrange = c(0, 0.8)),
  list(xrange = c(0, 1), yrange = c(0, 5)),
  list(xrange = c(0,5), yrange = c(0, 0.4))   
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
    select(-c(draws))
  list(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
       method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk'),
       data = temp)
}

boxplot_slice_sampler <- function(list_Of_df, curve_num, ksFilter) {
  
  ksFilterCutOff <- 0.10
  
  curveNames <- c('.2 norm(0,0.5) + .8 norm(6,2)','.2 norm(0,0.5) + .8 norm(20,1)','gamma(2.5,1)','truncated t(5)','t(5)','beta(0.2,0.8)','norm(20,5)')

  stepping_out <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'stepping_out')
  if(ksFilter == TRUE) {
    stepping_out <- stepping_out[[1]]$data %>% filter(ksTest >= ksFilterCutOff) %>%    mutate(parameters = paste('x =',x,'w =',w))
  } else {
    stepping_out <- stepping_out[[1]]$data %>% mutate(parameters = paste('x =',x,'w =',w))
  }
  
  # stepping out
  stepping_out_x <- stepping_out %>% 
    ggplot(aes(y = as.factor(x), x = SampPSec)) +
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  stepping_out_w <- stepping_out %>% 
    ggplot(aes(y = as.factor(w), x = SampPSec)) + 
    geom_boxplot() + 
    labs(
      x = 'Samples Per CPU Second',
      y = 'w'
    )+
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  stepping_out_param <- stepping_out %>% 
    ggplot(aes(y = as.factor(parameters), x = SampPSec)) + 
    geom_boxplot() + 
    labs(
      x = 'Samples Per CPU Second',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  stepping_out_samplesThin  <- stepping_out %>% 
    ggplot(aes(y = as.factor(parameters), x = samplesThin)) + 
    geom_boxplot() + 
    labs(
      x = 'Number Of Samples After Thinning',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) + 
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
    
  # latent
  
  latent <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'latent')
  if(ksFilter == TRUE) {
    latent <- latent[[1]]$data %>% filter(ksTest >= ksFilterCutOff) %>% mutate(parameters = paste('x = ',x,'s =',s,'rate =',rate))
  } else {
    latent <- latent[[1]]$data %>% mutate(parameters = paste('x = ',x,'s =',s,'rate =',rate))
  }

  latent_x <- latent %>% 
    ggplot(aes(y = as.factor(x), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  latent_s <- latent %>% 
    ggplot(aes(y = as.factor(s), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 's'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  latent_rate <- latent %>% 
    ggplot(aes(y = as.factor(rate), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'rate'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  latent_param <- latent %>% 
    ggplot(aes(x = SampPSec, y = as.factor(parameters))) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  latent_samplesThin  <- latent %>% 
    ggplot(aes(y = as.factor(parameters), x = samplesThin)) + 
    geom_boxplot() + 
    labs(
      x = 'Number Of Samples After Thinning',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) + 
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  # gess
  gess <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'gess')
  if(ksFilter == TRUE) {
    gess <- gess[[1]]$data %>% filter(ksTest >= ksFilterCutOff) %>%     mutate(parameters = paste('x =',x,'mu =',mu,'sigma =',sigma,'df =',df))
  } else {
    gess <- gess[[1]]$data %>% mutate(parameters = paste('x =',x,'mu =',mu,'sigma =',sigma,'df =',df))
  }

  gess_x <- gess %>% 
    ggplot(aes(y = as.factor(x), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  gess_mu <- gess %>% 
    ggplot(aes(y = as.factor(mu), x = SampPSec)) +
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'mu'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  gess_sigma <- gess %>% 
    ggplot(aes(y = as.factor(sigma), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'sigma'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  gess_df <- gess %>% 
    ggplot(aes(y = as.factor(df), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'df'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  gess_param <- gess %>% 
    ggplot(aes(x = SampPSec, y = as.factor(parameters))) +
    geom_boxplot() + 
    labs(
      x = 'Samples Per CPU Second',
      y = 'parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  gess_samplesThin  <- gess %>% 
    ggplot(aes(y = as.factor(parameters), x = samplesThin)) + 
    geom_boxplot() + 
    labs(
      x = 'Number Of Samples After Thinning',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) + 
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) +
    NULL
  
  
  # transform
  transform <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'transform')
  if(ksFilter == TRUE) {
    transform <- transform[[1]]$data %>% filter(ksTest >= ksFilterCutOff)
  } else {
    transform <- transform[[1]]$data
  }

  transform$inv_cdf_chr <- as.character(transform$inv_cdf) %>% str_extract('\\bq.*')
  transform <- transform %>% mutate(parameters = paste('x =',x,'inv_cdf =',inv_cdf_chr))
  
  transform_x <- transform %>% 
    ggplot(aes(y = as.factor(x), x = SampPSec)) + 
    geom_boxplot()  +
    labs(
      x = 'Samples Per CPU Second',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  transform_inv_cdf <- transform %>% 
    ggplot(aes(y = as.factor(inv_cdf_chr), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'Inverse CDF'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  transform_kld <- transform %>%
    mutate(parameters = case_when(
      method == 'transform' ~ stringr::str_remove_all(stringr::str_extract(parameters,'q.*'),'u,|=|mean|min|max|sd|q|shape1|shape2|shape|rate|df| '),
      TRUE ~ parameters
    )) %>% 
    mutate(parameters = case_when(
      grepl('*fit*',parameters) ~ 'LA',
      TRUE ~ parameters
    )) %>%
    ggplot(aes(y = SampPSec, x = KLD, color = parameters)) + 
    geom_point() +
    scale_y_continuous(labels = scales::number_format(scale = .001, suffix = 'K', accuracy = .1)) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      title = element_text(size=20),
      legend.title = element_text(size=20),
      legend.text = element_text(size=20),
      legend.position="right"
    ) +
    labs(
    title = curveNames[curve_num]
    ) +
    NULL
  
  
  transform_param <- transform %>% 
    ggplot(aes(y = as.factor(parameters), x = SampPSec)) +
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  transform_samplesThin  <- transform %>% 
    ggplot(aes(y = as.factor(parameters), x = samplesThin)) + 
    geom_boxplot() + 
    labs(
      x = 'Number Of Samples After Thinning',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) + 
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  # random walk
  rand_walk <- list.filter(list_Of_df, curve == curve_num) %>% 
    list.filter(method == 'rand_walk')
  if(ksFilter == TRUE) {
    rand_walk <- rand_walk[[1]]$data %>% filter(ksTest >= ksFilterCutOff) %>% mutate(parameters = paste('x =',x,'c =',c))
  } else {
    rand_walk <- rand_walk[[1]]$data %>%  mutate(parameters = paste('x =',x,'c =',c))
  }
  
  rand_walk_x <- rand_walk %>% 
    ggplot(aes(y = as.factor(x), x = SampPSec)) +
    geom_boxplot()  +
    labs(
      x = 'Samples Per CPU Second',
      y = 'x'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  rand_walk_c <- rand_walk %>% 
    ggplot(aes(y = as.factor(c), x = SampPSec)) +
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'c'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  
  rand_walk_param <- rand_walk %>% 
    ggplot(aes(y = as.factor(parameters), x = SampPSec)) + 
    geom_boxplot() +
    labs(
      x = 'Samples Per CPU Second',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  rand_walk_samplesThin  <- rand_walk %>% 
    ggplot(aes(y = as.factor(parameters), x = samplesThin)) + 
    geom_boxplot() + 
    labs(
      x = 'Number Of Samples After Thinning',
      y = 'Parameters'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) + 
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
fulldata <- rbind(stepping_out %>% select(method, parameters, SampPSec, ESS),
                  latent %>% select(method, parameters, SampPSec, ESS),
                  gess %>% select(method, parameters, SampPSec, ESS),
                  transform %>% select(method, parameters, SampPSec, ESS),
                  rand_walk %>% select(method, parameters, SampPSec, ESS))  
  
best_param <- fulldata[fulldata$method != 'transform',] %>%
  group_by(parameters, method) %>%
  summarise(avg = mean(SampPSec)) %>%
  arrange(desc(avg)) %>% 
  ungroup() %>% 
  group_by(method) %>% 
  slice_head(n = 1) %>% 
  select(parameters)

transform_params <- fulldata[fulldata$method == 'transform',] %>% 
  group_by(parameters) %>% 
  summarise(avg = mean(SampPSec)) %>%
  arrange(desc(avg)) 

if(max(grepl('*fit*', transform_params$parameters)) == 1) {
  laplaceApproximation <- transform_params[grepl('*fit*', transform_params$parameters),'parameters']
}

transform_params <- transform_params[round(quantile(1:nrow(transform_params), probs = c(0,0.33,0.66,1))),'parameters']
  
if(exists('laplaceApproximation')) {
  fulldata <- fulldata[fulldata$parameters %in% c(best_param$parameters, transform_params$parameters, laplaceApproximation$parameters),] %>% 
    mutate(parameters = ifelse(method != 'transform',method,parameters))
} else {
  fulldata <- fulldata[fulldata$parameters %in% c(best_param$parameters, transform_params$parameters),] %>% 
    mutate(parameters = ifelse(method != 'transform',method,parameters))
}

  # total comp
  sampler_comp <- fulldata %>% 
    mutate(parameters = case_when(
      method == 'transform' ~ stringr::str_remove_all(stringr::str_extract(parameters,'q.*'),'u,|=|mean|min|max|sd|q|shape1|shape2|shape|rate|df| '),
      TRUE ~ parameters
    )) %>% 
    mutate(parameters = case_when(
      grepl('*fit*',parameters) ~ 'Laplace Approximation',
      TRUE ~ parameters
    )) %>% 
    mutate(parameters = case_when(
      parameters == 'transform' ~ 'Transform',
      parameters == 'stepping_out' ~ 'Neal Stepping Out',
      parameters == 'rand_walk' ~ 'Random Walk',
      parameters == 'latent' ~ 'Latent',
      parameters == 'gess' ~ 'GESS',
      TRUE ~ parameters
    )) %>% 
    ggplot(aes(y = as.factor(parameters), x = SampPSec, fill = as.factor(method))) +
    geom_boxplot()  +
    labs(
      x = 'Samples Per CPU Second',
      y = 'Parameters',
      fill = 'Method'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  
  sampler_comp_samplesThin <- fulldata %>% 
    ggplot(aes(y = as.factor(parameters), x = ESS, fill = as.factor(method))) +
    geom_boxplot()  +
    labs(
      x = 'Number Of Samples After Thinning',
      y = 'Parameters',
      fill = 'Method'
    ) +
    scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) + 
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) + 
    NULL
  

  return(list(
    stepping_out_x = stepping_out_x, stepping_out_w = stepping_out_w, stepping_out_param = stepping_out_param, stepping_out_samplesThin = stepping_out_samplesThin,
    # stepping_out_plot = ((stepping_out_x + stepping_out_w)/(stepping_out_param)),
    latent_x = latent_x, latent_s = latent_s, latent_rate = latent_rate, latent_param = latent_param, latent_samplesThin = latent_samplesThin,
    # latent_plot = ((latent_x + latent_s)/(latent_rate)/(latent_param)),
    gess_x = gess_x, gess_mu = gess_mu, gess_sigma = gess_sigma, gess_df = gess_df, gess_param = gess_param, gess_samplesThin = gess_samplesThin,
    # gess_plot = ((gess_x + gess_mu)/(gess_sigma + gess_df)/(gess_param)),
    transform_x = transform_x, transform_inv_cdf = transform_inv_cdf, transform_param = transform_param, transform_samplesThin = transform_samplesThin, transform_kld = transform_kld,
    # transform_plot = ((transform_x) / (transform_inv_cdf) / (transform_kld) / (transform_param)),
    rand_walk_x = rand_walk_x, rand_walk_c = rand_walk_c, rand_walk_param = rand_walk_param, rand_walk_samplesThin = rand_walk_samplesThin,
    sampler_comp = sampler_comp,
    # rand_walk_plot = ((rand_walk_x + rand_walk_c)/(rand_walk_param)),
    curve_plot = (sampler_comp / sampler_comp_samplesThin)
  ))
}


extract_metrics <- function(path) {
  temp <- readRDS(path)
  temp <- temp %>% mutate(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
                          method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk')) %>% 
    select(method, curve, ESS, SampPSec, ksTest)
}



extract_draws <- function(path) {
  temp <- readRDS(path)
  temp <- temp %>% 
    mutate(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
           method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk'))
  if(unique(temp$method) == 'stepping_out') {
    temp <- temp %>% mutate(params = paste0(x,'|',w))
  } else if (unique(temp$method) == 'latent') {
    temp <- temp %>% mutate(params = paste0(x,'|',s,'|',rate))
  } else if (unique(temp$method) == 'gess') {
    temp <- temp %>% mutate(params = paste0(x,'|',mu,'|',sigma,'|',df))
  } else if (unique(temp$method) == 'transform') {
    temp <- temp %>% mutate(params = paste0(x,'|',str_extract(deparse1(inv_cdf), 'q.*')))
  } else if (unique(temp$method) == 'rand_walk') {
    temp <- temp %>% mutate(params = paste0(x,'|',c))
  }

  temp <- temp %>% 
    select(curve, method, params, draws, ksTest)
  
  list(curve = as.numeric(gsub(".*?([0-9]+).*", "\\1", path)),
       method = str_extract(path,'stepping_out|latent|gess|transform|rand_walk'),
       data = temp)
}


curve_slice_sampler <- function(list_Of_df, curve, ksFilter = FALSE) {
  
  ksFilterCutOff <- 0.10
  
  # curves and ranges
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
    list(xrange = c(-4, 15), yrange = c(0, 0.17)),
    list(xrange = c(-3, 25), yrange = c(0, 0.32)),
    list(xrange = c(0, 10), yrange = c(0, 0.31)),
    list(xrange = c(0, 10), yrange = c(0, 0.8)),
    list(xrange = c(-10, 10), yrange = c(0, 0.42)),
    list(xrange = c(0, 1), yrange = c(0, 5)),
    list(xrange = c(0, 40), yrange = c(0, 0.08)),
    list(xrange = c(0,5), yrange = c(0, 0.4))   
  )
  
  lf <- curves[[curve]]
  
  par(mfrow = c(2,3))
  # stepping out
  stepping_out <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'stepping_out') 
  if(ksFilter == TRUE) {
    stepping_out <- stepping_out[[1]]$data %>% filter(ksTest >= ksFilterCutOff)
  } else {
    stepping_out <- stepping_out[[1]]$data
  }
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2
  )
  title(main = 'Stepping Out')
  lapply(stepping_out$draws, function(x) {
    lines(density(x, adjust = 2), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2,
    add = TRUE
  )
  
  # latent 
  latent <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'latent')
  if(ksFilter == TRUE) {
    latent <- latent[[1]]$data %>% filter(ksTest >= ksFilterCutOff)
  } else {
    latent <- latent[[1]]$data
  }
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2
  )
  title(main = 'Latent')
  lapply(latent$draws, function(x) {
    lines(density(x, adjust = 2), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2,
    add = TRUE
  )
  
  # gess
  gess <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'gess')
  if(ksFilter == TRUE) {
    gess <- gess[[1]]$data %>% filter(ksTest >= ksFilterCutOff)
  } else {
    gess <- gess[[1]]$data
  }
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2
  )
  title(main = 'GESS')
  lapply(gess$draws, function(x) {
    lines(density(x, adjust = 2), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2,
    add = TRUE
  )
  
  # transform
  transform <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'transform')
  if(ksFilter == TRUE) {
    transform <- transform[[1]]$data %>% filter(ksTest >= ksFilterCutOff)
  } else {
    transform <- transform[[1]]$data
  }
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2
  )
  title(main = 'Transform')
  lapply(transform$draws, function(x) {
    lines(density(x, adjust = 2), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2,
    add = TRUE
  )
  
  # rand_walk
  rand_walk <- list.filter(list_Of_df, curve == curve) %>% 
    list.filter(method == 'rand_walk')
  if(ksFilter == TRUE) {
    rand_walk <- rand_walk[[1]]$data %>% filter(ksTest >= ksFilterCutOff)
  } else {
    rand_walk <- rand_walk[[1]]$data
  }
  
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2
  )
  title(main = 'Random Walk')
  lapply(rand_walk$draws, function(x) {
    lines(density(x, adjust = 2), col = adjustcolor('black'))
  })
  curve(
    fexp(f = lf, x = x),
    col = 'red',
    xlim = range[[curve]]$xrange,
    ylim = range[[curve]]$yrange,
    xlab = '',
    ylab = '',
    lwd = 2,
    add = TRUE
  )
}



facet_plot <- function(list) {
  lf <- curves[[i]]
  xlim <- range[[i]]$xrange
  ylim <- range[[i]]$yrange
  splitList <- split(list[[1]]$data, f = list[[1]]$data$params)
  par(mfrow = c(ceiling(length(splitList)/3),3), mar = c(1,1,1,1))
  for(i in 1:length(splitList)) {
    tempDf <- splitList[[i]]
    targetRate <- round(mean(tempDf$ksTest >= 0.05),2) * 100
    ksTestColor <- ifelse(tempDf$ksTest >= 0.05, 'dodgerblue','firebrick')
    curve(
      fexp(f = lf, x = x),
      col = 'black',
      xlim = xlim,
      ylim = ylim,
      xlab = '',
      ylab = '',
      lwd = 2,
      xaxt = 'n',
      yaxt = 'n',
      main = paste0(tempDf$params %>% unique(),' ', targetRate, '%'),
    )
    for (i in 1:nrow(tempDf)) {
      lines(density(tempDf$draws[[i]]), col = ksTestColor[i])
    }
  }
}

