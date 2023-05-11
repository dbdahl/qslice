# loading in libraries
library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)
source("../sam_comparison_functions.R")


reps <- 1

fexp <- function(f, x) exp(f(x))
samples <- rep(500000,reps)
sampleSize <- 5000


# this function finds the second derivative of function f at point x
second_derivative <- function( x, h = 1e-5, f ) {
  
  num <- f(x + h) - 2*f(x) + f(x - h)
  denom <- h^2
  
  num/denom
}

# this function provides a lapace approximation for a curve lf
laplace_approx <- function(lf, h = 1e-5, init, inflate = TRUE){
  fit <- optim(init, lf, control=list(fnscale=-1), method = 'BFGS')
  if(inflate){
    sd <- sqrt(-solve(second_derivative(fit$par, h = h, f = lf))) * 3
  } else {
    sd <- sqrt(-solve(second_derivative(fit$par, h = h, f = lf)))
  }
  return(
      c(log_pdf = function(x) {dnorm(x, mean = fit$par, sd = sd, log = TRUE)},
      inv_cdf = function(u, lower.tail = TRUE) {qnorm(u, mean = fit$par, sd = sd, lower.tail = lower.tail)},
      mu = fit$par,
      sigma = sd)
  )
}

# this function provides a cauchy approximation for a curve lf
cauchy_approx <- function(lf, h = 1e-5, init, scaleInflation = 1) {
  temp <- laplace_approx(lf, h = h, init = init, inflate = FALSE)
  return(
    c(log_pdf = function(x) {dcauchy(x, location = temp$mu, scale = temp$sigma)},
      inv_cdf = function(u) {qcauchy(u, location = temp$mu, scale = temp$sigma)},
      loc = temp$mu,
      scale = scaleInflation*temp$sigma)
  )
}




# helpful functions
# truncated cauchy
dcauchy_trunc <- function(x, loc, sc, lb=-Inf, ub=Inf, log=FALSE) {
  stopifnot(lb < ub)
  
  denom <- pcauchy(ub, loc=loc, scale=sc) - pcauchy(lb, loc=loc, scale=sc)
  
  if (log) {
    out <- dcauchy(x, loc=loc, sc=sc, log=TRUE) - log(denom)
    out[which(x <= lb)] = -Inf
    out[which(x >= ub)] = -Inf
  } else {
    out <- dcauchy(x, loc=loc, sc=sc) / denom
    out[which(x <= lb)] = 0.0
    out[which(x >= ub)] = 0.0
  }
  
  out
}

# give draws and will return scale and location for a cauchy
fit_trunc_Cauchy <- function(y, lb=-Inf, ub=Inf) {
  
  llik <- function(params, y, lb, ub) {
    loc <- params[1]
    sc <- params[2]
    
    stopifnot(all(y < ub))
    stopifnot(all(y > lb))
    
    if (sc > 0.0) {
      llik <- dcauchy_trunc(y, loc=loc, sc=sc, lb=lb, ub=ub, log=TRUE)
      out <- sum(llik)
    } else {
      out <- -Inf
    }
    
    out
  }
  
  fit <- optim(par = c(mean(y), sd(y)), fn=llik, control=list(fnscale=-1),
               y = y, lb = lb, ub = ub)
  
  list(loc = fit$par[1], sc = fit$par[2], lb = lb, ub = ub, fit = fit)
}

# this function will return the log pdf and inv cdf for a cauchy
pseudo_Cauchy <- function(loc, sc, lb=-Inf, ub=Inf) {
  
  plb <- pcauchy(lb, loc=loc, sc=sc)
  pub <- pcauchy(ub, loc=loc, sc=sc)
  normc <- pub - plb
  
  pseudo_log_pdf <- function(x) dcauchy(x, loc=loc, sc=sc, log=TRUE) - log(normc)
  pseudo_inv_cdf <- function(u) qcauchy(plb + u*normc) * sc + loc
  
  list(pseudo_log_pdf = pseudo_log_pdf,
       pseudo_inv_cdf = pseudo_inv_cdf,
       loc = loc, sc = sc)
}

# this function returns a Laplace t approximation with log pdf and inv cdf
lapproxt <- function(lf, init, sc_adj = 1.0, lb = -Inf, ub = Inf, ...) {
  
  fit <- optim(par = init, fn = lf, control = list(fnscale = -1), method = 'BFGS')
  loc <- fit$par
  hessian <- second_derivative( x = loc, h = 1e-5, f = lf )
  hessian_f <- hessian * exp(fit$value) # hessian of original f
  sc <- sc_adj / sqrt(-hessian_f)
  out <- pseudo_Cauchy(loc = loc, sc = sc, lb = lb, ub = ub)
  out[["fit"]] <- fit
  
  out
}
