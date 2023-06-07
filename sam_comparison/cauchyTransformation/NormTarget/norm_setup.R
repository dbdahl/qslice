## Set up for the cauch transform 
## author: Sam Johnnson

# setwd("~/cucumber/sam_comparison")
source('../../sam_comparison_functions.R')
library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)
library(rootSolve)

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
pseudo_Cauchy <- function(loc, sc, lb=-Inf, ub=Inf, name = NULL) {
  
  plb <- pcauchy(lb, loc=loc, sc=sc)
  pub <- pcauchy(ub, loc=loc, sc=sc)
  normc <- pub - plb
  
  pseudo_log_pdf <- function(x) dcauchy(x, loc=loc, sc=sc, log=TRUE) - log(normc)
  pseudo_inv_cdf <- function(u) qcauchy(plb + u*normc) * sc + loc
  
  list(pseudo_log_pdf = pseudo_log_pdf,
       pseudo_inv_cdf = pseudo_inv_cdf,
       t = paste0("Cauchy(loc=",round(loc,2), ", sc=",round(sc,2),"), ", name),
       loc = loc, sc = sc)
}


# returns a psuedo target given a function
lapproxt <- function(f, init, sc_adj = 1.0, lb = -Inf, ub = Inf, minlb = -1e6, maxub = 1e6, ...) {
  
  tempUb <- ifelse(is.infinite(ub), maxub, ub)
  tempLb <- ifelse(is.infinite(lb), minlb, lb)
  loc <- optimize(f, interval = c(tempLb,tempUb), maximum = TRUE)$maximum
  while (is.infinite(loc) | any(abs(c(tempLb - loc, tempUb - loc)) < 1)) {
    if(abs(tempLb - loc) < 1) tempLb <- tempLb * .80
    if(abs(tempUb - loc) < 1) tempUb <- tempUb * .80
    loc <- optimize(f, interval = c(tempLb,tempUb), maximum = TRUE)$maximum
  }
  # getting the second derivative  
  hessian <- numDeriv::hessian(func = f, x = loc, method = 'Richardson')[1,1]
  sc <- sc_adj / sqrt(-hessian)
  out <- pseudo_Cauchy(loc = loc, sc = sc, lb = lb, ub = ub)
  # out[["fit"]] <- fit # this is needed if the maximium is found using optim
  out <- append(out, list(hessian = hessian))
  
  out
}

## optimizes the area under the curve to find loc and scale

pseudo_Cauchy_list <-  function(loc, sc, lb=-Inf, ub=Inf) {
  
  plb = pcauchy(lb, loc=loc, sc=sc)
  pub = pcauchy(ub, loc=loc, sc=sc)
  normc = pub - plb
  
  list(d = function(x) {dcauchy(x, loc=loc, sc=sc)},
       ld = function(x) {dcauchy(x, loc=loc, sc=sc, log=TRUE)},
       dld = function(x) {-2*(x-loc)/sc^2 / (1 + ((x-loc)/sc)^2)}, # derivative of log density
       q = function(u) {qcauchy(plb + u*normc)*sc + loc},
       p = function(x) {(pcauchy(x, loc=loc, sc=sc) - plb) / normc},
       t = paste0("Cauchy(", loc, ", ", sc, ")"))
}

# samples = rnorm(10e3)
# samples = rgamma(10e3, 2.5, 1)

opt_Cauchy_auc_data <-  function(samples, lb=-Inf, ub=Inf) {
  
  get_auc = function(pars, samples, lb, ub) {
    loc = pars[1]
    sc = pars[2]
    
    pseu = pseudo_Cauchy_list(loc=loc, sc=sc, lb=lb, ub=ub)
    
    qq = pseu$p(samples)
    
    nbins = 30
    (bins = seq(0.0, 1.0, len=nbins+1))
    (tab = tabulate( as.numeric(cut(qq, breaks = bins)), nbins=nbins))
    
    tab[c(1,nbins)] = 1.2 * tab[c(1,nbins)] # penalty for smiling...
    
    (tab_norm = tab / max(tab) / nbins)
    
    (auc = sum(tab_norm))
    auc
  }
  
  temp <- optim(c(0.0, 1.0), get_auc, control = list(fnscale=-1), samples=samples,
                lb = lb, ub = ub)
  
  # temp$par[1] <- ifelse(temp$par[1] < lb, lb, temp$par[1])
  # temp$par[1] <- ifelse(temp$par[2] > ub, ub, temp$par[1])
  
  list(loc = temp$par[1], sc = temp$par[2], lb = lb, ub = ub)
  
}

# opt_Cauchy_auc_data(samples, lb=-Inf, ub=Inf)
# opt_Cauchy_auc_data(samples, lb=0, ub=Inf)

opt_Cauchy_auc <-  function(target) {
  
  get_auc = function(pars, target) {
    loc = pars[1]
    sc = pars[2]
    
    pseu = pseudo_Cauchy_list(loc=loc, sc=sc, lb=target$lb, ub=target$ub)
    
    lh = function(u, targ, pseu) {
      (targ$ld( pseu$q(u) ) - pseu$ld( pseu$q(u) ))
    }
    
    h = function(u, targ, pseu) {
      lh(u, targ, pseu) |> exp()
    }
    
    dlh = function(x, targ, pseu) {
      targ$dld(x) - pseu$dld(x)
    }
    
    (extrema_x = uniroot.all(dlh, interval=c(pseu$q(0.05), pseu$q(0.95)),
                             n = 5e3,
                             targ=target, pseu=pseu))
    
    (extrema_u = pseu$p(extrema_x))
    
    (extrema_h = h(extrema_u, targ=target, pseu=pseu))
    
    m_indx = which.max(extrema_h)
    n_extrema = length(extrema_x)
    (m = extrema_h[m_indx])
    
    (tmp = c(0.0, extrema_u, 1.0))
    (mids_u = tmp[-length(tmp)] + diff(tmp)/2)
    (mids_x = pseu$q(mids_u))
    (dlh_at_mids = dlh(mids_x, targ=target, pseu=pseu))
    
    (local_max = ((dlh_at_mids[m_indx] > 0) && (dlh_at_mids[m_indx + 1] < 0)))
    
    auc = integrate(function(u, t, p){h(u, targ=t, pseu=p) / m}, lower=0.0, upper=1.0, t=target, p=pseu)
    stopifnot(auc$message == "OK")
    
    pen_dip = 0.0 * ifelse(length(extrema_x) == 1, 0.0, max(extrema_h) - min(extrema_h))
    # pen_dip = ifelse(length(extrema_x) == 1, 0.0, 1.0)
    
    pen_loc_min = 1.0 * !local_max
    
    auc$value - pen_dip - pen_loc_min
  }
  
  optim(c(0.0, 1.2), get_auc, control = list(fnscale=-1), target = target)
  
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

scales <- c(0.75,0.875,1,1.25,1.35,1.5)

makePseudo <- function(cauchy_scale) {
  temp <- lapproxt(f = \(x) exp(lf(x)), init = 1, sc_adj = cauchy_scale, minlb = -3, maxub = 3)
  # temp <- pseudo_Cauchy(loc = temp$loc, sc = temp$sc)
  list(
    d = temp$pseudo_log_pdf,
    q = temp$pseudo_inv_cdf,
    t = paste0("Cauchy(loc=",round(temp$loc,2), ", sc=",round(temp$sc,2),")")
  )
}

### Specialty Scenarios
## dr heiner optimization
optimizedCauchy <- pseudo_Cauchy(loc = 0.05328, sc = 1.2702, name = 'Optim Samps')
# standard loc 0, sc 1
standardCauchy <- pseudo_Cauchy(loc = 0, sc = 1)
# fit using Auto
# getting burn in draws to fit pseudo target
burnin_metrics = stepping_out_time_eval(
  samples = min(5000, max(samples) * 0.1),
  x_0 = x,
  lf_func = lf,
  w_value = median(w),
  max_value = Inf,
  log_value = TRUE
)

# fitting the Cauchy
psuedoFit <- fit_trunc_Cauchy(unlist(burnin_metrics$Draws))
autoCauchy <- pseudo_Cauchy(loc = psuedoFit$loc,sc = psuedoFit$sc, name = 'Auto')
# fit using Optim Samples
psuedoFit <- opt_Cauchy_auc_data(unlist(burnin_metrics$Draws))
optimSamplesCauchy <- pseudo_Cauchy(loc = psuedoFit$loc,sc = psuedoFit$sc, name = 'Optim Samps')
# fit using optim
truth = list(d = function(x) {dnorm(x)},
             ld = function(x) {dnorm(x, log=TRUE)},
             # dld = function(x) {-x},
             q = function(u) {qnorm(u)},
             lb = -Inf, ub = Inf,
             t = "normal(0,1)")
# truth$dld = Deriv::Deriv(truth$ld, x="x")
truth$dld = function(x) numDeriv::grad(truth$ld, x=x)
optimCauchy <- opt_Cauchy_auc(truth)
optimCauchy <- pseudo_Cauchy(loc = optimCauchy$par[1], sc = optimCauchy$par[2], name = 'Optim')

## transform tuning parameters ##
log_pdf <- lapply(scales, FUN = \(par) makePseudo(par)) %>% 
  lapply(FUN = \(x) x$d) %>% 
  unlist()


log_pdf <- append(log_pdf, list(optimizedCauchy$pseudo_log_pdf, standardCauchy$pseudo_log_pdf,
                  autoCauchy$pseudo_log_pdf, optimSamplesCauchy$pseudo_log_pdf, optimCauchy$pseudo_log_pdf))

inv_cdf <- lapply(scales, FUN = \(par) makePseudo(par)) %>% 
  lapply(FUN = \(x) x$q) %>% 
  unlist()

inv_cdf <- append(inv_cdf, list(optimizedCauchy$pseudo_inv_cdf, standardCauchy$pseudo_inv_cdf,
                  autoCauchy$pseudo_inv_cdf, optimSamplesCauchy$pseudo_inv_cdf, optimCauchy$pseudo_inv_cdf))

t <- lapply(scales, FUN = \(par) makePseudo(par)) %>% 
  lapply(FUN = \(x) x$t) %>% 
  unlist()

t <- c(t, list(optimizedCauchy$t, standardCauchy$t, 
       autoCauchy$t, optimSamplesCauchy$t, optimCauchy$t))

## random walk tuning parameters ##
c <- c(0.25,0.5,1,1.5,2,2.5,3,3.5,4)
