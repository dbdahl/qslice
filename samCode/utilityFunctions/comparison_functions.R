# rm(list=ls())
library(progress)
library(magrittr)
library(tidyverse)
library(utils)
# library(cucumber)

############ Stepping Out Eval #############


# function to evaluate stepping out procedure
stepping_out_time_eval <- function(samples = 200,
                                   x_0 = 0.5,
                                   lf_func,
                                   w_value = .5,
                                   max_value = Inf,
                                   log_value = TRUE,
                                   prog_bar = TRUE) {
  # browser()
  counter <- 0
  draws <- numeric(samples)
  draws[1] <- x_0
  if(prog_bar) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  time <- system.time({
    print(paste0(x_0,"|",w_value))
    for ( i in seq.int(2,length(draws)) ) {
      out <- slice_sampler_stepping_out(x = draws[i-1], target=lf_func, w=w_value, max=max_value, log=log_value)
      draws[i] <- out$x
      counter <- counter + out$nEvaluations
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS,
                         userTime = time['user.self'], sysTime = time['sys.self'], elapsedTime = time['elapsed'])
  temp_tbl$Draws <- list(draws)
  temp_tbl
}




############ Latent Eval ############


# creating a function to evaluate the latent slice sampler
latent_time_eval <- function(samples = 50000,
                             x_0,
                             s_0,
                             lf_func,
                             rate_value,
                             log_value = TRUE,
                             prog_bar = TRUE) {
  # browser()
  counter <- 0
  draws <- latents <- numeric(samples)
  draws[1] <- x_0
  latents[1] <- s_0
  if(prog_bar) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  time <- system.time({
    print(paste0(x_0,"|",s_0,"|",rate_value))
    for ( i in seq.int(2,length(draws)) ) {
    # i <- 2
    # while(coda::effectiveSize(draws) < samples) {
      out <- slice_sampler_latent(draws[i-1], latents[i-1], target = lf_func, rate=rate_value)
      draws[i] <- out$x
      latents[i] <- out$s
      counter <- counter + out$nEvaluations
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
      # i <- i + 1
    }
  })
  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS,
                         userTime = time['user.self'], sysTime = time['sys.self'], elapsedTime = time['elapsed'])
  temp_tbl$Draws <- list(draws)
  temp_tbl
}


################### Elliptical Eval ##################


## Samples with elliptical slice sampler (Murray 2010)
elliptical_time_eval <- function(samples = 50000,
                             lf_func,
                             x_0 = 2,
                             mu_value = 2,
                             sigma_value = 5,
                             log_value = TRUE,
                             prog_bar = TRUE) {
  counter <- 0
  draws <- matrix(nrow = samples, ncol = 2)
  draws[1,] <- x_0
  if(prog_bar) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  time <- system.time({
    print(paste0(x_0,"|",mu_value,"|",sigma_value))
    for ( i in seq.int(2,nrow(draws)) ) {
      out <- slice_sampler_elliptical(x = matrix(draws[i-1,], ncol = 2),
                                                mu = mu_value,
                                                sigma = sigma_value,
                                                target = lf_func,
                                                log = log_value)
      draws[i,] <- out$x
      counter <- counter + out$nEvaluations
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS,
                         userTime = time['user.self'], sysTime = time['sys.self'], elapsedTime = time['elapsed'])
  temp_tbl$Draws <- list(draws)
  temp_tbl
}

################ GESS ##############

## Samples with generalized elliptical slice sampler (Nishihara 2014)
gess_time_eval <- function(samples = 50000,
                           lf_func,
                           x_0 = 2,
                           mu_value = 2,
                           sigma_value = 5,
                           df_value = 4,
                           upr_support = Inf,
                           lwr_support = -Inf,
                           log_value = TRUE,
                           prog_bar = TRUE){
  # browser()
  counter <- 0
  draws <- numeric(samples)
  draws[1] <- x_0
  if(prog_bar) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  time <- system.time({
    print(paste0(x_0,"|",mu_value,"|",sigma_value,"|",df_value))
    for ( i in seq.int(2,length(draws)) ) {
      out <- slice_sampler_generalized_elliptical(x = draws[i-1], mu = mu_value, sigma = sigma_value, target = lf_func, df = df_value, log_value)
      # checking to see if x is in the support if not redraw
      while(out$x >= upr_support | out$x <= lwr_support) {
          print("outside of support")
          print(out$x)
          print(i)
          out <- slice_sampler_generalized_elliptical(x = draws[i-1], mu = mu_value, sigma = sigma_value, target = lf_func, df = df_value, log_value)
        }
      draws[i] <- out$x
      counter <- counter + out$nEvaluations
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS,
                         userTime = time['user.self'], sysTime = time['sys.self'], elapsedTime = time['elapsed'])
  temp_tbl$Draws <- list(draws)
  temp_tbl
}


############## Transform Eval ################

# function to evaluate transform procedure
transform_time_eval <- function(samples = 200,
                                   x_0 = 0.5,
                                   lf_func,
                                   pseudo_pdf_log,
                                   pseudo_cdf_inv,
                                   log_value = TRUE,
                                   prog_bar = TRUE) {
  # browser()
  counter <- 0
  draws <- numeric(samples)
  Tdraws <- numeric(samples)
  draws[1] <- x_0
  Tdraws[1] <- 0
  if(prog_bar) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  time <- system.time({
    print(paste0(x_0,"|",deparse1(pseudo_cdf_inv)))
    for ( i in seq.int(2,length(draws)) ) {
      out <- slice_sampler_transform(x = draws[i-1], target=lf_func, pseudo_log_pdf = pseudo_pdf_log, pseudo_inv_cdf = pseudo_cdf_inv, log=log_value)
      draws[i] <- out$x
      Tdraws[i] <- out$u
      counter <- counter + out$nEvaluations
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS,
                         userTime = time['user.self'], sysTime = time['sys.self'], elapsedTime = time['elapsed'])
  temp_tbl$Draws <- list(draws)
  temp_tbl$TDraws <- list(Tdraws)
  temp_tbl
}


############## Random Walk ################

# function to evaluate Random Walk
random_walk_time_eval <- function(samples = 200,
                         x_0 = 0.5,
                         c = 2,
                         support = c(-Inf,Inf),
                         lf_func,
                         prog_bar = TRUE) {
  ## setup to count number of evaluations
  nEvaluations <- 0
  f <- function(x) { nEvaluations <<- nEvaluations + 1; lf_func(x) }
  #MCMC setup
  # counter <- 0
  draws <- numeric(samples) #initialize a vector to save the accepted values of the parameter
  draws[1] <- x_0
  int.x <- x_0 #starting value of parameter
  n.accept <- 0 #how many times do we accept the proposed value?
  if(prog_bar) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  #MCMC algorithm (metropolis random walk)
  time <- system.time({
    print(paste0(int.x,"|",c))
    for(i in 2:samples){
      # proposed draw
      x.dot <- rnorm(1, int.x, c)
      if(x.dot >= support[1] && x.dot <= support[2]){

          logr <- f(x.dot) - f(int.x)

          u <- runif(1, 0, 1)
          if(log(u) < logr){
            int.x <- x.dot
            n.accept <- n.accept+1
          }
      }
      draws[i] <- int.x
      # counter <- counter + 1
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
    }
  })
  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_random_walk <- time['user.self']
  temp_tbl <- data.frame(nEval = samples, EffSamp = ESS, acceptanceRate = n.accept/samples,
                         userTime = time['user.self'], sysTime = time['sys.self'], elapsedTime = time['elapsed'])
  
  temp_tbl$Draws <- list(draws)
  temp_tbl
}


