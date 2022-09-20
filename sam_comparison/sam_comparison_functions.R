rm(list=ls())
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
    for ( i in seq.int(2,length(draws)) ) {
      out <- slice_sampler_stepping_out(x = draws[i-1], target=lf_func, w=w_value, max=max_value, log=log_value)
      # checking to see if x is in the support and redrawing
      # while(out$x >= upr_support | out$x <= -Inf) {
      #   print("outside of support")
      #   print(out$x)
      #   print(i)
      #   out <- slice_sampler_generalized_elliptical(x = draws[i-1], mu = mu_value, sigma = sigma_value, target = lf_func, df = df_value, log_value)
      # }
      draws[i] <- out$x
      counter <- counter + out$nEvaluations
      if(prog_bar) {utils::setTxtProgressBar(pb, i)}
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_step_out <- time['user.self']
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS, Time = time_step_out)
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
                             log_value = TRUE) {
  # browser()
  counter <- 0
  draws <- latents <- numeric(samples)
  draws[1] <- x_0
  latents[1] <- s_0
  pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  time <- system.time({
    for ( i in seq.int(2,length(draws)) ) {
      out <- slice_sampler_latent(draws[i-1], latents[i-1], target = lf_func, rate=rate_value)
      # checking to see if x is in the support if not redrawing
      # while(out$x >= upr_support | out$x <= -Inf) {
      #   print("outside of support")
      #   print(out$x)
      #   print(i)
      #   out <- slice_sampler_generalized_elliptical(x = draws[i-1], mu = mu_value, sigma = sigma_value, target = lf_func, df = df_value, log_value)
      # }
      draws[i] <- out$x
      latents[i] <- out$s
      counter <- counter + out$nEvaluations
      utils::setTxtProgressBar(pb, i)
    }
  })
  # draws_latent <- draws
  # counter_latent <- counter
  # time_latent <- coda::effectiveSize(draws) / time['user.self']

  # counter
  # plot(density(draws))
  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_latent <- time['user.self']
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS, Time = time_latent)
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
                             pb = TRUE) {
  counter <- 0
  draws <- matrix(nrow = samples, ncol = 2)
  draws[1,] <- x_0
  if(pb) {pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)}
  time <- system.time({
    for ( i in seq.int(2,nrow(draws)) ) {
      out <- slice_sampler_elliptical(x = matrix(draws[i-1,], ncol = 2),
                                                mu = mu_value,
                                                sigma = sigma_value,
                                                target = lf_func,
                                                log = log_value)
      # checking to see if x is in the support and if not redrawing
      # while(out$x >= upr_support | out$x <= -Inf) {
      #   print("outside of support")
      #   print(out$x)
      #   print(i)
      #   out <- slice_sampler_generalized_elliptical(x = draws[i-1], mu = mu_value, sigma = sigma_value, target = lf_func, df = df_value, log_value)
      # }
      draws[i,] <- out$x
      counter <- counter + out$nEvaluations
      if(pb) {utils::setTxtProgressBar(pb, i)}
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_elliptical <- time['user.self']
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS, Time = time_elliptical)
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
                           log_value = TRUE){
  # browser()
  counter <- 0
  draws <- numeric(samples)
  draws[1] <- x_0
  pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  time <- system.time({
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
      utils::setTxtProgressBar(pb, i)
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_elliptical <- time['user.self']
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS, Time = time_elliptical)
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
                                   log_value = TRUE) {
  # browser()
  counter <- 0
  draws <- numeric(samples)
  draws[1] <- x_0
  pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  time <- system.time({
    for ( i in seq.int(2,length(draws)) ) {
      out <- slice_sampler_transform(x = draws[i-1], target=lf_func, pseudo_log_pdf = pseudo_pdf_log, pseudo_inv_cdf = pseudo_cdf_inv, log=log_value)
      # checking to see if x is in the support and redrawing
      # while(out$x >= upr_support | out$x <= -Inf) {
      #   print("outside of support")
      #   print(out$x)
      #   print(i)
      #   out <- slice_sampler_generalized_elliptical(x = draws[i-1], mu = mu_value, sigma = sigma_value, target = lf_func, df = df_value, log_value)
      # }
      draws[i] <- out$x
      counter <- counter + out$nEvaluations
      utils::setTxtProgressBar(pb, i)
    }
  })

  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_transform <- time['user.self']
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS, Time = time_transform)
  temp_tbl$Draws <- list(draws)
  temp_tbl
}


############## Random Walk ################

# function to evaluate Random Walk
random_walk_time_eval <- function(samples = 200,
                         x_0 = 0.5,
                         c = 2,
                         support = c(-Inf,Inf),
                         lf_func) {
  #MCMC setup
  counter <- 0
  draws <- numeric(samples) #initialize a vector to save the accepted values of the parameter
  int.x <- x[1] #starting value of parameter
  n.accept <- 0 #how many times do we accept the proposed value?
  pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  #MCMC algorithm (metropolis random walk)
  time <- system.time({
    for(i in 1:samples){
      # proposed draw
      x.dot <- rnorm(1, int.x, c)
      if(x.dot >= support[1] && x.dot <= support[2]){

          logr <- lf_func(x.dot) - lf_func(int.x)

          u <- runif(1, 0, 1)
          if(log(u) < logr){
            int.x <- x.dot
            n.accept <- n.accept+1
          }
      }
      draws[i] <- int.x
      counter <- counter + 1
      utils::setTxtProgressBar(pb, i)
    }
  })
  ESS <- ifelse(coda::effectiveSize(coda::as.mcmc(draws)) > samples,
                samples,
                coda::effectiveSize(coda::as.mcmc(draws)))
  time_random_walk <- time['user.self']
  temp_tbl <- data.frame(nEval = counter, EffSamp = ESS, Time = time_random_walk)
  temp_tbl$Draws <- list(draws)
  temp_tbl
}


# ## Samples with slice_sampler_transform
# ## Pseudo prior (select one)
# pseudoLogPDF <- function(x) dnorm(x, mean=4.0, sd=15.0, log=TRUE)
# pseudoInvCDF <- function(u) qnorm(u, mean=4.0, sd=15.0)
#
# pseudoLogPDF <- function(x) dt((x-ellip_ctr)/ellip_sc, df=ellip_degf, log=TRUE) - log(ellip_sc)
# pseudoInvCDF <- function(u) ellip_sc * qt(u, df=ellip_degf) + ellip_ctr
#
# pseudoLogPDF <- function(x) dgamma(x, shape=2.5, log=TRUE)
# pseudoInvCDF <- function(u) qgamma(u, shape=2.5)
#
# pseudoLogPDF <- function(x) dexp(x, log=TRUE)
# pseudoInvCDF <- function(u) qexp(u)
#
# pseudoLogPDF <- function(x) if ( x < 0 ) -Inf else dt(x, df=1.0, log=TRUE) + log(2.0)
# pseudoInvCDF <- function(u) qt((u + 1.0)/2.0, df=1) # half Cauchy
#
# pseudoLogPDF <- function(x) dbeta(x, shape1=0.5, shape2=0.5, log=TRUE)
# pseudoInvCDF <- function(u) qbeta(u, shape1=0.5, shape2=0.5)
#
#
# counter <- 0
# draws <- numeric(50000)
# time <- system.time({
#   for ( i in seq.int(2,length(draws)) ) {
#     out <- slice_sampler_transform(x = draws[i-1], pseduo_log_pdf = pseduoLogPDF, pseduo_inv_cdf = pseduoInvCDF, log_density = lf) #x, log_density, pseudo_log_pdf, pseudo_inv_cdf
#     draws[i] <- out$x
#     counter <- counter + out$nEvaluations
#   }
# })
# draws_transform <- draws
# counter_transform <- counter
# plot(density(draws))
# time_generalized_elliptical <- coda::effectiveSize(draws) / time['user.self']
#
#
#
####### The code below is not currently tested / functional ##############
#
#
#
# # Setup for elliptical slice sampler
# ellip_ctr <- 0.0
# ellip_sc <- 1.0
# ellip_degf <- 10.0
#
# ellip_invTform <- function(z) z # identity transformation
# ellip_logJacobian <- function(z) 0.0
#
# ellip_invTform <- function(z) exp(z) # log transformation
# ellip_logJacobian <- function(z) z
#
# ellip_invTform <- function(z) 1.0 / (1.0 + exp(-z)) # logit transformation
# ellip_logJacobian <- function(z) {
#   enz <- exp(-z)
#   -z - 2.0*log(1.0 + enz)
# }
#
# ellip_logtarget <- function(z) lf( ellip_invTform(z) ) + ellip_logJacobian(z)
#
#
# ## Pseudo prior (select one)
# pseudoLogPDF <- function(x) dnorm(x, mean=4.0, sd=15.0, log=TRUE)
# pseudoInvCDF <- function(u) qnorm(u, mean=4.0, sd=15.0)
#
# pseudoLogPDF <- function(x) dt((x-ellip_ctr)/ellip_sc, df=ellip_degf, log=TRUE) - log(ellip_sc)
# pseudoInvCDF <- function(u) ellip_sc * qt(u, df=ellip_degf) + ellip_ctr
#
# pseudoLogPDF <- function(x) dgamma(x, shape=2.5, log=TRUE)
# pseudoInvCDF <- function(u) qgamma(u, shape=2.5)
#
# pseudoLogPDF <- function(x) dexp(x, log=TRUE)
# pseudoInvCDF <- function(u) qexp(u)
#
# pseudoLogPDF <- function(x) if ( x < 0 ) -Inf else dt(x, df=1.0, log=TRUE) + log(2.0)
# pseudoInvCDF <- function(u) qt((u + 1.0)/2.0, df=1) # half Cauchy
#
# pseudoLogPDF <- function(x) dbeta(x, shape1=0.5, shape2=0.5, log=TRUE)
# pseudoInvCDF <- function(u) qbeta(u, shape1=0.5, shape2=0.5)
#
#
#
#
# ## Samples with proposed method
# counter <- 0
# lf_compos <- function(u) {
#   x <- pseudoInvCDF(u)
#   lf(x) - pseudoLogPDF(x)
# }
# draws1 <- numeric(50000)
# time1 <- system.time({
#   u <- 0.5
#   lfx <- lf_compos(u)
#   if ( lfx == -Inf ) stop("Oops, the starting value is not in support.")
#   current <- list(x=u, lfx=lfx)
#   for ( i in seq_along(draws1) ) {
#     current <- slice_sampler_shrinkage(current, lf_compos)
#     draws1[i] <- pseudoInvCDF(current$x)
#   }
# })
# counter1 <- counter
#
# counter <- 0
# draws10 <- numeric(50000)
# draws10[1] <- 0.5
# time10 <- system.time({
#   for ( i in seq_along(draws1)[-1] ) {
#     draws10[i] <- slice_sampler_transform(draws10[i-1], lf, pseudoLogPDF, pseudoInvCDF)
#   }
# })
# counter10 <- counter
#
#
# ## Samples with latent slice method (Li & Walker, 2020)
# counter <- 0
# draws2 <- numeric(50000)
# time2 <- system.time({
#   x <- 0.1
#   s <- 1.0
#   lambda <- 0.1
#   lfx <- lf(x)
#   if ( lfx == -Inf ) stop("Oops, the starting value is not in support.")
#   current <- list(x=x, s=s, lfx=lfx)
#   for ( i in seq_along(draws1) ) {
#     current <- slice_latent(current, lf, lambda)
#     draws2[i] <- current$x
#   }
# })
# counter2 <- counter
#
#
# ## Samples with elliptical slice sampler
# draws3 <- numeric(50000)
# time3 <- system.time({
#   z <- 0.1
#   lfz <- ellip_logtarget(z)
#   if ( lfz == -Inf ) stop("Oops, the starting value is not in support.")
#   current <- list(x=z, lfx=lfz)
#   for ( i in seq_along(draws3) ) {
#     current <- slice_ellipse_generalized(current, ellip_logtarget, ellip_ctr, ellip_sc, ellip_degf)
#     draws3[i] <- ellip_invTform(current$x)
#   }
# })
# counter3 <- counter
#
# ## version of elliptical sampler with no transform (considerably faster than identity transform)
# # draws2 <- numeric(50000)
# # time2 <- system.time({
# #   z <- 0.1
# #   lfz <- lf(z)
# #   if ( lfz == -Inf ) stop("Oops, the starting value is not in support.")
# #   current <- list(x=z, lfx=lfz)
# #   for ( i in seq_along(draws2) ) {
# #     current <- slice_ellipse_generalized(current, lf, ellip_ctr, ellip_sc, ellip_degf)
# #     draws2[i] <- current$x
# #   }
# # })
# # counter2 <- counter
# # counter <- 0
#
#
#
#
# ## check
# # Stepping out
# hist(draws0[-c(1:10000)], freq=FALSE, breaks=50)
# lines(density(draws0), col="red")
# curve(f, -3, 25, add=TRUE, col="black", n=1000)
# (B0 <- coda::effectiveSize(draws0))   # Effective number of independent samples
# B0 / time0['user.self']             # Effective number of independent samples per CPU second
# counter0                           # Number of function evaluations
#
# x_test = 0.7
#
# integrate(f, -Inf, x_test)
# mean(draws0 < x_test)
# (err0 <- abs( integrate(f, -Inf, x_test)$value - mean(draws0 < x_test) ))
#
# (acf0 <- acf(draws0)$acf[2,1,1])
#
#
# # Transform
# hist(draws1[-c(1:10000)], freq=FALSE, breaks=50)
# lines(density(draws1), col="red")
# curve(f, -3, 25, add=TRUE, col="black", n=1000)
# curve(exp(pseudoLogPDF(x)), -3, 25, add=TRUE, col="blue", lty=2)
# (B1 <- coda::effectiveSize(draws1))  # Effective number of independent samples
# B1 / time1['user.self']             # Effective number of independent samples per CPU second
# counter1                            # Number of function evaluations
#
# integrate(f, -Inf, x_test)
# mean(draws1 < x_test)
# (err1 <- abs( integrate(f, -Inf, x_test)$value - mean(draws1 < x_test) ))
#
# (acf1 <- acf(draws1)$acf[2,1,1])
#
#
# # Latent
# hist(draws2[-c(1:10000)], freq=FALSE, breaks=50)
# lines(density(draws2), col="red")
# curve(f, -3, 25, add=TRUE, col="black", n=1000)
# (B2 <- coda::effectiveSize(draws2))  # Effective number of independent samples
# B2 / time2['user.self']             # Effective number of independent samples per CPU second
# counter2                            # Number of function evaluations
#
# integrate(f, -Inf, x_test)
# mean(draws2 < x_test)
# (err2 <- abs( integrate(f, -Inf, x_test)$value - mean(draws2 < x_test) ))
#
# (acf2 <- acf(draws2)$acf[2,1,1])
#
#
# # Elliptical
# hist(draws3[-c(1:10000)], freq=FALSE, breaks=50)
# lines(density(draws3), col="red")
# curve(f, -3, 25, add=TRUE, col="black", n=1000)
# curve( dt((x-ellip_ctr)/ellip_sc, df=ellip_degf, log=TRUE) - log(ellip_sc), -3, 25, add=TRUE, col="blue", lty=2) # only makes sense if no transformation
# (B3 <- coda::effectiveSize(draws3))  # Effective number of independent samples
# B3 / time3['user.self']             # Effective number of independent samples per CPU second
# counter3                            # Number of function evaluations
#
# integrate(f, -Inf, x_test)
# mean(draws3 < x_test)
# (err3 <- abs( integrate(f, -Inf, x_test)$value - mean(draws3 < x_test) ))
#
# (acf3 <- acf(draws3)$acf[2,1,1])
#
#
#
#
#
# ## compare
# err0 # stepping out
# err1 # transform
# err2 # latent
# err3 # t-elliptical
#
# acf0
# acf1
# acf2
# acf3
#
# (B0 <- coda::effectiveSize(draws0))   # stepping out
# (B1 <- coda::effectiveSize(draws1))   # transform
# (B10 <- coda::effectiveSize(draws10))   # transform
# (B2 <- coda::effectiveSize(draws2))   # latent
# (B3 <- coda::effectiveSize(draws3))   # t-elliptical
#
# B0 / time0['user.self'] # stepping out
# B1 / time1['user.self'] # transform
# B10 / time10['user.self'] # transform
# B2 / time2['user.self'] # latent
# B3 / time3['user.self'] # t-elliptical
#
# counter0
# counter1
# counter2
# counter3

