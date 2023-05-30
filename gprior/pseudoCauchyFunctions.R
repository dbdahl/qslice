## function to fit a psuedo cauchy
## Author: Matthew Heiner and Sam Johnson

library(rootSolve)
library(Deriv)
library(numDeriv)

Nburnin <- 1000

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

# get samples using stepping out procedure
burnin <- function(x, lf, samples, w) {
  draws <- numeric(samples)
  draws[1] <- x
  for( i in 2:samples ) {
    if(sum(is.na(draws)) != 0) browser()
    draws[i] <- cucumber::slice_sampler_stepping_out(draws[i-1], target = lf, w = w, log = TRUE, max = Inf)$x
  }
  draws
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

# given a function returns the second derivative at that point
second_derivative <- function( x, h = 1e-7, f ) {
  
  num <- f(x + h) - 2*f(x) + f(x - h)
  # num <- f(x + 2*h) - 2*f(x + h) + f(x)
  # num <- f(x) - 2*f(x - h) + f(x - 2*h)
  denom <- h^2
  
  num/denom
}

# returns a psuedo target given a function
lapproxt <- function(f, init, sc_adj = 1.0, lb = -Inf, ub = Inf, ...) {
  
  fit <- optim(par = init, fn = f, control = list(fnscale = -1), method = 'BFGS')
  loc <- fit$par
  hessian <- second_derivative( x = loc, h = 1e-7, f = f )
  # hessian_f <- hessian * exp(fit$par) # hessian of original f
  # if(hessian_f > 0 ) browser()
  # sc <- sc_adj / sqrt(-hessian_f)
  # if(hessian > 0) browser()
  sc <- sc_adj / sqrt(-hessian)
  out <- pseudo_Cauchy(loc = loc, sc = sc, lb = lb, ub = ub)
  out[["fit"]] <- fit
  out <- append(out, list(hessian = hessian))
  
  out
}

## optimizes the area under the curve to find loc and scale

pseudo_Cauchy_list = function(loc, sc, lb=-Inf, ub=Inf) {
  
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

opt_Cauchy_auc_data = function(samples, lb=-Inf, ub=Inf) {
  
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

opt_Cauchy_auc = function(target) {
  
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


## Standard normal
truth = list(d = function(x) {dnorm(x)},
             ld = function(x) {dnorm(x, log=TRUE)},
             # dld = function(x) {-x},
             q = function(u) {qnorm(u)},
             lb = -Inf, ub = Inf,
             t = "normal(0,1)")
truth$dld = Deriv::Deriv(truth$ld, x="x")
truth$dld = function(x) numDeriv::grad(truth$ld, x=x)

opt_Cauchy_auc(truth)

## seeing how much area under the curve
auc_diagnostic = function(samples_u, nbins = 30) {
  
  (bins = seq(0.0, 1.0, len=nbins+1))
  (tab = tabulate( as.numeric(cut(samples_u, breaks = bins)), nbins=nbins))
  
  tab[c(1,nbins)] = 1.2 * tab[c(1,nbins)] # penalty for smiling...
  
  (tab_norm = tab / max(tab) / nbins)
  
  (auc = sum(tab_norm))
  auc
}
