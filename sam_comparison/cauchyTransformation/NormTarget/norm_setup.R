## Set up for the cauch transform 
## author: Sam Johnnson

# setwd("~/cucumber/sam_comparison")
source('../../sam_comparison_functions.R')
library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)

reps <- 10

fexp <- function(f, x) exp(f(x))
samples <- rep(5e4,reps)
saveInd <- TRUE

# functions used to create a psuedo target
# the pdf for a truncated cauchy
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

# given draws this function returns a fitted cauchy
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

# given location and scale returns a pdf and inverse cdf on the log scale
pseudo_Cauchy <- function(loc, sc, lb=-Inf, ub=Inf) {
  
  plb <- pcauchy(lb, loc=loc, sc=sc)
  pub <- pcauchy(ub, loc=loc, sc=sc)
  normc <- pub - plb
  
  pseudo_log_pdf <- function(x) dcauchy(x, loc=loc, sc=sc, log=TRUE) - log(normc)
  pseudo_inv_cdf <- function(u) qcauchy(plb + u*normc) * sc + loc
  
  list(pseudo_log_pdf = pseudo_log_pdf,
       pseudo_inv_cdf = pseudo_inv_cdf,
       t = paste0("Cauchy(loc=",round(loc,2), ", sc=",round(sc,2),"),Auto"),
       loc = loc, sc = sc)
}


second_derivative <- function( x, h = 1e-5, f ) {
  
  num <- f(x + h) - 2*f(x) + f(x - h)
  denom <- h^2
  
  num/denom
}

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


# 1
################### TEST ########################
## dnorm(x,0,1) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  dnorm(x,0,1,log = TRUE)
}


# making a grid to calculate KL divergence
grid <- seq(from = qnorm((1-.99999)/2, mean = 0, sd = 1),
            to = qnorm((1-.99999)/2, mean = 0, sd = 1, lower.tail = FALSE),
            length.out = 5000)

py <- exp(lf(grid))

xlim_range <- c(-4, 15)
ylim_range <- c(0, 0.42)

#### Tuning Parameters ####
## starting point ##
x <- c(0)#c(0, 1, 5)

## stepping out metrics to input ##
w <- c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.05, 0.1, 0.5, 1, 2, 5)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(0)#c(1,2,3,4.5,6,7)
sigma <- c(0.25, 0.5,1, 2, 3)#c(2,3,4,5,6,8)
df <- c(3, 5)#c(1,4,16,16^2,16^4)

scales <- c(0.5,0.75,0.875,1,1.25,1.35,1.5,2)

makePseudo <- function(cauchy_scale) {
  temp <- lapproxt(lf = lf, init = 1, sc_adj = cauchy_scale)
  # temp <- pseudo_Cauchy(loc = temp$loc, sc = temp$sc)
  list(
    d = temp$pseudo_log_pdf,
    q = temp$pseudo_inv_cdf,
    t = paste0("Cauchy(loc=",round(temp$loc,2), ", sc=",round(temp$sc,2),")")
  )
}

## transform tuning parameters ##
log_pdf <- lapply(scales, FUN = \(par) makePseudo(par)) %>% 
  lapply(FUN = \(x) x$d) %>% 
  unlist()

inv_cdf <- lapply(scales, FUN = \(par) makePseudo(par)) %>% 
  lapply(FUN = \(x) x$q) %>% 
  unlist()

t <- lapply(scales, FUN = \(par) makePseudo(par)) %>% 
  lapply(FUN = \(x) x$t) %>% 
  unlist()

## random walk tuning parameters ##
c <- c(0.25,0.5,1,1.5,2,2.5,3,3.5,4)
