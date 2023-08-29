# source("fn_AUC_hist.R")
# source("fn_skew_hist.R")
# source("fn_waterArea_hist.R")
# source("fn_waterArea_int.R")
# source("fn_expectedSliceWidth.R")


## version with AUC
utility_shrinkslice <- function(h = NULL, x = NULL, y = NULL, u = NULL,
                                type = "samples", # supplied with u. Alternatively, type = "samples_kde", type = "function" (supplied with function h) or "grid" (supplied with x, y)
                                coeffs = c(1.0, 1.0),
                                nbins = 30,
                                plot = FALSE, 
                                use_meanSliceWidth = TRUE,
                                tol_int = 1.0e-3) {

  if (type %in% c("function", "samples_kde")) {

    ## the supplied function here is the transformed target with support on (0, 1)

    if (type == "samples_kde") {
      h = beta_kde(u)
    }

    tmp <- water_area_int(h, interval = c(0.0, 1.0), plot = plot, eps = 1.0e-3, tol_int = tol_int)
    auc <- tmp$AUC / tmp$totalArea
    wtr <- tmp$totalWaterArea / tmp$totalArea
    
    if (use_meanSliceWidth) {
      msw <- meanSliceWidth_int(h, tol = tol_int)
    }

  } else if (type %in% c("grid", "samples")) {
    
    if (type == "samples") {
      
      if (is.null(x)) {
        x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins) # trouble if it doesn't reach far enough into tails
      } else {
        stopifnot(length(x) == nbins)
      }
      
      bins <- c(0.0, x[-c(1, nbins)], 1.0)
      y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins = nbins )
      if(sum(y) == 0) {
        if (use_meanSliceWidth) {
          out <- c( util = 0, auc = 0, water_area = 0, msw = 0 )
        } else {
          out <- c( util = 0, auc = 0, water_area = 0 )
        }
        return(out)
      } 
    }

    auc <- auc(x = x, y = y)
    wtr <- water_area(x = x, y = y, plot = plot)
    
    if (use_meanSliceWidth) {
      msw <- meanSliceWidth_grid(x = x, y = y)
    }
    
  } # currently no support for histogram method (unless supplied as x and y)

  if (use_meanSliceWidth) {
    out <- c( util = coeffs %*% c(msw, -wtr), auc = auc, water_area = wtr, msw = msw )
  } else {
    out <- c( util = coeffs %*% c(auc, -wtr), auc = auc, water_area = wtr )
  }
  
  out
}



# x = seq(1.0e-6, 1.0 - 1.0e-6, length=30)
# # x = runif(50, 0, 1)
# 
# use_msw <- TRUE

# dev.off()
# par(mfrow=c(1,2))
# y = dbeta(x, 0.5, 0.5); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 0.5, 0.5), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = dbeta(x, 1.0, 0.8); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 1.0, 0.8), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = dbeta(x, 1.0, 1.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 1.0, 1.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = dbeta(x, 2.0, 2.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 2.0, 2.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = dbeta(x, 1.0, 2.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 1.0, 2.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = dbeta(x, 20.0, 20.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 20.0, 20.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = dbeta(x, 3.0, 8.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) dbeta(x, 3.0, 8.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = 0.9*dbeta(x, 3.0, 8.0) + 0.1*dbeta(x, 15, 1.5); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) 0.9*dbeta(x, 3.0, 8.0) + 0.1*dbeta(x, 15, 1.5), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = 0.5*(x >= 0.5)*dbeta(x, 1.0, 2.0) + 0.5*(x < 0.5)*dbeta(x, 2.0, 1.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) 0.5*(x >= 0.5)*dbeta(x, 1.0, 2.0) + 0.5*(x < 0.5)*dbeta(x, 2.0, 1.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = 0.5*(x < 0.5)*(1.0 - 2.0*x) + 0.5*(x >= 0.5)*2.0*(x - 0.5); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) 0.5*(x < 0.5)*(1.0 - 2.0*x) + 0.5*(x >= 0.5)*2.0*(x - 0.5), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = 0.4*dbeta(x, 30.0, 10.0) + 0.6*dbeta(x, 10.0, 30.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) 0.4*dbeta(x, 30.0, 10.0) + 0.6*dbeta(x, 10.0, 30.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = 0.4*dbeta(x, 30.0, 3.0) + 0.6*dbeta(x, 3.0, 30.0) + 0.2; utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) 0.4*dbeta(x, 30.0, 3.0) + 0.6*dbeta(x, 3.0, 30.0) + 0.2, type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = 2.0*dbeta(x, 30.0, 10.0) + 5.0*dbeta(x, 3, 3) + 1.7*dbeta(x, 10.0, 30.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw); utility_shrinkslice(h = \(x) 2.0*dbeta(x, 30.0, 10.0) + 5.0*dbeta(x, 3, 3) + 1.7*dbeta(x, 10.0, 30.0), type = "function", plot = TRUE, use_meanSliceWidth = use_msw)
# y = runif(length(x), min = 0.6, max = 1.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = use_msw)
# dev.off()
# 

# pseudo_Cauchy_list <- function(loc, sc, lb = -Inf, ub = Inf, log_p = FALSE) {
#   
#   plb <- pcauchy(lb, loc=loc, sc=sc)
#   pub <- pcauchy(ub, loc=loc, sc=sc)
#   normc <- pub - plb
#   
#   list(d = function(x) {dcauchy(x, loc=loc, sc=sc)},
#        ld = function(x) {dcauchy(x, loc=loc, sc=sc, log=TRUE)},
#        dld = function(x) {-2*(x-loc)/sc^2 / (1 + ((x-loc)/sc)^2)},
#        q = function(u, log.p = FALSE) {qcauchy(plb + u*normc, log.p = log.p)*sc + loc},
#        p = function(x) {(pcauchy(x, loc=loc, sc=sc) - plb) / normc},
#        t = paste0("Cauchy(", loc, ", ", sc, ")"))
# }

pseudo_t_list <- function(loc, sc, degf, lb = -Inf, ub = Inf, log_p = FALSE, name = NULL) {
  
  if(!is.null(name)) {
    t <- paste0("t(loc = ", round(loc,2), ", sc = ", round(sc,2), ", degf = ", round(degf), "), ", name)
  } else {
    t <- paste0("t(loc = ", round(loc,2), ", sc = ", round(sc,2), ", degf = ", round(degf), ")")
  }
  
  plb <- pt((lb - loc)/sc, df = degf)
  pub <- pt((ub - loc)/sc, df = degf)
  normc <- pub - plb
  
  logsc <- log(sc)
  
  list(d = function(x) {dt((x - loc)/sc, df = degf) / sc},
       ld = function(x) {dt((x - loc)/sc, df = degf, log=TRUE) - logsc},
       # dld = function(x) {-2*(x-loc)/sc^2 / (1 + ((x-loc)/sc)^2)},
       q = function(u, log.p = FALSE) {qt(plb + u*normc, log.p = log.p, df = degf)*sc + loc},
       p = function(x) {(pt((x - loc)/sc, df = degf) - plb) / normc},
       t = t,
       loc = loc, sc = sc, degf = degf)
}


util_pseu <- function(pseu, target = NULL, samples = NULL, 
                      type = "samples", 
                      x = NULL, 
                      bins = NULL, nbins = NULL, 
                      coeffs = c(1.0, 1.0), 
                      plot = FALSE, 
                      use_meanSliceWidth = TRUE,
                      tol_int = 1.0e-3) {
  
  if (type == "function") {
    
    h <- function(x) exp( target$ld( pseu$q(x) ) - pseu$ld( pseu$q(x) ) )
    x <- NULL
    y <- NULL
    u <- NULL
    
  } else {
    
    h <- NULL
    
    if (type == "grid") {
      
      xx <- pseu$q(x)
      ly <- target$ld( xx ) - pseu$ld( xx )
      y <- exp(ly - max(ly))
      u <- NULL
      
    } else if (type == "samples") {
      
      u <- pseu$p(samples) # if type == "samples_kde"
      # y <- tabulate( as.numeric(cut(qq, breaks = bins)), nbins = nbins )
      y <- NULL
      
    } else if (type == "samples_kde") {
      
      u <- pseu$p(samples) # if type == "samples_kde"
      # y <- tabulate( as.numeric(cut(qq, breaks = bins)), nbins = nbins )
      x <- NULL
      y <- NULL
      
    }
    
  }
  
  utility_shrinkslice(h = h, x = x, y = y, u = u,
                      type = type,
                      coeffs = coeffs, 
                      nbins = nbins, 
                      plot = plot, 
                      use_meanSliceWidth = use_meanSliceWidth,
                      tol_int = tol_int)
}


opt_t <- function(target = NULL, samples = NULL, 
                  type = "samples", # one of "samples", "grid", "function"
                  degf = c(1, 5, 20), 
                  lb = -Inf, ub = Inf, 
                  nbins = 100,
                  coeffs = c(1.0, 1.0),
                  tol_opt = 1.0e-6, tol_int = 1.0e-3,
                  plot = TRUE,
                  verbose = FALSE,
                  use_meanSliceWidth = TRUE) {
  
  if (type %in% c("function", "samples_kde")) {
    x <- NULL
    bins <- NULL
    nbins <- NULL
  } else if (type %in% c("grid", "samples")) {
    x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins) # trouble if it doesn't reach far enough into tails
    bins <- c(0.0, x[-c(1, nbins)], 1.0)
  }
  
  if (is.null(samples)) {
    inits <- c(loc = 0.5, sc = 2.0)
  } else {
    inits <- c(loc = mean(samples), sc = sd(samples))
  }
  
  get_util <- function(pars, type, x, target, samples, degf, lb, ub, bins, nbins, 
                       coeffs, verbose, use_meanSliceWidth) {
    
    loc <- pars[1]
    sc <- pars[2]
    
    if (sc <= 0.0) {
      
      out <- -10.0
      
    } else {
      
      pseu <- pseudo_t_list(loc = loc, sc = sc, degf = degf,
                            lb = lb, ub = ub)
      
      if (verbose) cat("trying", pseu$t, "\n")
      
      out <- util_pseu(pseu = pseu, target = target, samples = samples, 
                       type = type, x = x, 
                       bins = bins, nbins = nbins,
                       coeffs = coeffs, plot = FALSE,
                       use_meanSliceWidth = use_meanSliceWidth,
                       tol_int = tol_int)
      
      if (verbose) {
        print(out)
        cat("\n") 
      }
      
    }
    
    out["util"]
  }
  
  
  opt <- list()
  
  for (i in 1:length(degf)) {
    opt[[i]] <- optim(inits, get_util, control = list(fnscale=-1, reltol = tol_opt),
                      type = type,
                      x = x, target = target, samples = samples, 
                      degf = degf[i],
                      lb = lb, ub = ub,
                      bins = bins, nbins = nbins,
                      coeffs = coeffs, 
                      verbose = verbose,
                      use_meanSliceWidth = use_meanSliceWidth)
  }
  
  use_indx <- which.max(sapply(opt, function(obj) obj$value))
  
  pseu <- pseudo_t_list(loc = opt[[use_indx]]$par[1], 
                        sc = opt[[use_indx]]$par[2], 
                        degf = degf[use_indx],
                        lb = lb, ub = ub)
  
  util <- util_pseu(pseu = pseu, target = target, samples = samples, 
                    type = type, 
                    x = x, 
                    bins = bins, nbins = nbins,
                    coeffs = coeffs, 
                    plot = plot,
                    use_meanSliceWidth = use_meanSliceWidth,
                    tol_int = tol_int)
  
  list(pseu = pseu, util = util, opt = opt[[use_indx]], nbins = nbins, coeffs = coeffs, 
       tol_int = tol_int, tol_opt = tol_opt)
}



#### 
## Examples
####


# 
# coeffs <- c(1.0, 0.0) # opt auc loc = 0, sc = 1.0, degf = 20; 
# coeffs <- c(1.0, 0.0) # opt msw loc = 0, sc = 0.98, degf = 20
# coeffs <- c(1.0, 1.0) # opt auc - water loc = 0.00, sc = 1.0, degf = 20
# coeffs <- c(1.0, 1.0) # opt msw - water loc = 0.01, sc = 1.0, degf = 20
# coeffs <- c(1.0, 100) # opt auc - water loc = 0.004, sc = 1.015, degf = 20
# coeffs <- c(1.0, 100) # opt msw - water loc = 0.01, sc = 1.0, degf = 20
# 
# use_msw <- TRUE
# 
# samples <- rnorm(2e3)
# # (pseu <- opt_skewt(samples = samples, nbins = 100))
# # (pseu <- opt_skewt(target = list(ld = function(x) dnorm(x, log = TRUE)), nbins = 200))
# (pseu <- opt_t(samples = samples, nbins = 30, coeffs = coeffs, plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw))
# (pseu <- opt_t(target = list(ld = function(x) dnorm(x, log = TRUE)), type = "grid", nbins = 100, coeffs = coeffs, plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw, tol_opt = 1e-3))
# (pseu <- opt_t(target = list(ld = function(x) dnorm(x, log = TRUE)), type = "function", coeffs = coeffs, plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw, tol_opt = 1e-3, tol_int = 1e-2))
# 
# curve(dnorm(x), from = -3, to = 3)
# curve(pseu$pseu$d(x), from = -3, to = 3, add = T, col = "red")
# 
# curve( exp( dnorm(pseu$pseu$q(x), log=TRUE) - pseu$pseu$ld(pseu$pseu$q(x)) ), from = 0, to = 1 )
# 
# uu <- pseu$pseu$p(rnorm(2e3))
# # uu <- sn::pst(rnorm(500e3), xi = 0, omega = 1.0, nu = 99, alpha = 0.0)
# # uu <- runif(500e3)
# utility_shrinkslice(u = uu, nbins = 30, plot = TRUE) # not implemented
# 
# # nbins <- 100
# # x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins)
# # bins <- c(0.0, x[-c(1, nbins)], 1.0)
# # xx <- pseu$pseu$q(x)
# # ly <- dnorm(xx, log = TRUE) - pseu$pseu$ld( xx )
# # y <- exp(ly - max(ly))
# # utility_shrinkslice(x = x, y = y, nbins = nbins, plot = TRUE)
# # plot(xx, dnorm(xx))
# # points(xx, pseu$pseu$d(xx), col="red")
# # plot(x,y)
# # dnorm(xx) / pseu$pseu$d(xx)
# # log(dnorm(xx) / pseu$pseu$d(xx))
# # dnorm(xx, log=T) - pseu$pseu$ld(xx)
# 
# 
# coeffs <- c(1.0, 0.0) # opt auc loc = 1.47, sc = 1.81, degf = 5; water 0.02; auc .877
# coeffs <- c(1.0, 0.0) # opt msw loc = 1.73, sc = 1.69, degf = 5; water 0.06
# coeffs <- c(1.0, 1.0) # opt msw - water loc = 1.35, sc = 1.87, degf = 5; water 0.006; auc .870
# coeffs <- c(1.0, 1.0) # opt msw - water loc = 1.59, sc = 1.84, degf = 5
# coeffs <- c(1.0, 100) # opt msw - water loc = 1.30, sc = 1.93, degf = 5; auc .857
# coeffs <- c(1.0, 100) # opt msw - water loc = 1.56, sc = 1.89, degf = 5
# 
# pseu <- list(pseu = pseudo_t_list(loc = 1.30, sc = 1.93, degf = 5, lb = 0))
# 
# use_msw <- TRUE
# 
# samples <- rgamma(5e3, 2.5, 1)
# # (pseu <- opt_skewt(samples = samples, nbins = 30))
# # (pseu <- opt_skewt(target = list(ld = function(x) dgamma(x, 2.5, 1.0, log = TRUE)), nbins = 50, lb = 0.0))
# (pseu <- opt_t(samples = samples, nbins = 30, lb = 0.0, coeffs = coeffs, plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw))
# (pseu <- opt_t(target = list(ld = function(x) dgamma(x, 2.5, 1.0, log = TRUE)), type = "grid", nbins = 30, lb = 0.0, coeffs = coeffs, plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw, tol_opt = 1e-3))
# (pseu <- opt_t(target = list(ld = function(x) dgamma(x, 2.5, 1.0, log = TRUE)), type = "function", lb = 0.0, coeffs = coeffs, plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw, tol_opt = 1e-3, tol_int = 1e-2))
# 
# curve(dgamma(x, 2.5, 1), from = 0, to = 10)
# curve(pseu$pseu$d(x), from = 0, to = 10, add = T, col = "red")
# 
# curve( exp( dgamma(pseu$pseu$q(x), 2.5, 1, log=TRUE) - pseu$pseu$ld(pseu$pseu$q(x)) ), from = 0, to = 1 )
# curve( dgamma(pseu$pseu$q(x), 2.5, 1, log=FALSE) / pseu$pseu$d(pseu$pseu$q(x)), from = 0, to = 1 )
# 
# uu <- pseu$pseu$p(rgamma(50e3, 2.5, 1))
# utility_shrinkslice(u = uu, nbins = 100, plot = TRUE)
# 
# 
# 
# 
# coeffs <- c(1.0, 0.0) # opt auc loc = 0.34, sc = 0.41, degf = 1; water 0.03; auc .795
# coeffs <- c(1.0, 0.0) # opt msw loc = 0.41, sc = 0.38, degf = 1; water 0.08
# coeffs <- c(1.0, 1.0) # opt auc - water loc = 0.30, sc = 0.43, degf = 1; water 0.01; auc 0.789
# coeffs <- c(1.0, 1.0) # opt msw - water loc = 0.37, sc = 0.44, degf = 1; water 0.009; util 0.816
# coeffs <- c(1.0, 100) # opt auc - water loc = 0.27, sc = 0.45, degf = 1
# coeffs <- c(1.0, 100) # opt msw - water loc = 0.33, sc = 0.47, degf = 1
# 
# pseu <- list(pseu = pseudo_t_list(loc = 0.33, sc = 0.47, degf = 1, lb = 0))
# 
# use_msw = TRUE
# 
# source("fn_invgamma.R")
# 
# samples <- rinvgamma(50e3, 2, 1)
# # (pseu <- opt_skewt(samples = samples, nbins = 30))
# # (pseu <- opt_skewt(target = list(ld = function(x) dgamma(x, 2.5, 1.0, log = TRUE)), nbins = 50, lb = 0.0))
# (pseu <- opt_t(samples = samples, nbins = 30, lb = 0.0, coeffs = coeffs, degf = c(1, 5), plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw))
# (pseu <- opt_t(target = list(ld = function(x) dinvgamma(x, 2.0, 1.0, log = TRUE)), type = "grid", nbins = 30, lb = 0.0, coeffs = coeffs, degf = c(1), plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw, tol_opt = 1e-3))
# (pseu <- opt_t(target = list(ld = function(x) dinvgamma(x, 2.0, 1.0, log = TRUE)), type = "function", lb = 0.0, coeffs = coeffs, degf = c(1), plot = TRUE, verbose = TRUE, use_meanSliceWidth = use_msw, tol_opt = 1e-3, tol_int = 1e-2))
# 
# curve(dinvgamma(x, 2.0, 1.0), from = 0, to = 4)
# curve(pseu$pseu$d(x), from = 0, to = 4, add = T, col = "red")
# 
# curve( exp( dinvgamma(pseu$pseu$q(x), 2.0, 1.0, log=TRUE) - pseu$pseu$ld(pseu$pseu$q(x)) ), from = 0, to = 1 )
# curve( dinvgamma(pseu$pseu$q(x), 2.0, 1.0, log=FALSE) / pseu$pseu$d(pseu$pseu$q(x)), from = 0, to = 1 )
# 
# uu <- pseu$pseu$p(rinvgamma(50e3, 2.0, 1.0))
# utility_shrinkslice(u = uu, nbins = 100, plot = TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
# ### Sam's samples from g-prior
# load("test/approxSamples.rda")
# (pseu <- opt_t(samples = samples, nbins = 30, lb = 0.0))
# hist(samples, freq=FALSE)
# curve(pseu$pseu$d(x), from = 0, to = 200, add = T, col = "red")
# (pseu <- opt_t(samples = samples, nbins = 50, lb = 0.0, coeffs = c(1,1)))
# hist(samples, freq=FALSE)
# curve(pseu$pseu$d(x), from = 0, to = 200, add = T, col = "orange")
# hist(pseu$pseu$p(samples))
# utility_shrinkslice(u = pseu$pseu$p(samples), nbins = 30, plot = TRUE)
# 
