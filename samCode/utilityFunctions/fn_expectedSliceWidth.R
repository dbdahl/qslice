meanSliceWidth_int <- function(h, tol = 0.005, eps = 1e-3) {
  
  interval <- c(0.0,1.0)
  # if (is.infinite(h(interval[1])) || is.nan(h(interval[1]))) {
  #   interval[1] <- interval[1] + eps
  # } 
  # if (is.infinite(h(interval[2])) || is.nan(h(interval[2]))) {
  #   interval[2] <- interval[2] - eps
  # }
  # if (is.infinite(h(interval[1])) || is.nan(h(interval[1])) || is.infinite(h(interval[2])) || is.nan(h(interval[2]))) return(0)

  nc <- integrate(h, lower = interval[1], upper = interval[2])$value
  h_norm <- function(x) h(x) / nc
  
  h_fill <- function(x, hnorm_at_outer) { # x must be able to be a vector
    pmin(h_norm(x), hnorm_at_outer)
  }
  
  inner_int <- function(x_outer) {
    y <- h_norm(x_outer)
    sapply(y, function(z) { # inner_int must be vectorized
      integrate(h_fill, lower = interval[1], upper = interval[2], hnorm_at_outer = z, abs.tol = tol)$value
    })
  }
  
  integrate(inner_int, lower = interval[1], upper = interval[2], abs.tol = tol)$value
}


# source("fn_waterArea_int.R")
# 
# auc_w <- function(h) {
#   wai <- water_area_int(h, plot = TRUE)
#   c(wai$AUC / wai$totalArea, wai$totalWaterArea / wai$totalArea)
# }
# 



# x <- seq(0,1, by = 0.2)
# y <- runif(length(x))
# 
# plot(x,y)
# segments(x, 0, y1 = y)
# abline(h = y, lty = 2)
# 
# y
# (ord <- order(y, y + (1.0e-9)*runif(length(y), min = -1.0, max = 1.0)))
# (yord <- y[ord])
# (A_inYorder <- (length(yord):1)*yord + c(0, cumsum(yord[-length(yord)])))
# stopifnot(all.equal(yord[order(ord)], y))
# A_inYorder[order(ord)]


# decif <- function(boolvec) {
#   n <- length(boolvec)
#   out <- numeric(n)
#   out[1] <- n
#   saveup <- 0
#   for (i in 2:n) {
#     if (boolvec[i]) {
#       out[i] <- out[i-1] - 1 - saveup
#       saveup <- 0
#     } else {
#       out[i] <- out[i-1]
#       saveup <- saveup + 1
#     }
#   }
#   out  
# }
# 
# decif(c(T,T,F,F,T,T,F,T,F,F))
# 
# yord
# c(TRUE, diff(yord) > 0.0)
# decif(c(TRUE, diff(yord) > 0.0))


meanSliceWidth_grid <- function(x, y) {
  n <- length(x)
  
  stopifnot(length(y) == n)
  
  sumy <- sum(y)
  yn <- y / sum(y)
  
  ord <- order(y, y + (1.0e-9)*runif(n, min = -1.0, max = 1.0)) # to break ties... not necessary?
  yn_ord <- yn[ord]
  A_inYorder <- (n:1)*yn_ord + c(0, cumsum(yn_ord[-n]))
  
  mean(A_inYorder)
}


beta_kde <- function(samples, ...) {
  require("bde")
  
  dd <- bde::chen99Kernel(samples, ...)
  h <- function(x) density(dd, x)
  h
}

# h <- beta_kde(samples)
# curve(h, from = 0, to = 1)

 
# xx <- seq(1.0e-6, 1.0 - 1.0e-6, length = 1000)
# 
# h <- function(x) dbeta(x, 2, 1); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 1, 2); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 1, 1); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 2, 6); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dunif(x, min = 0.5, max = 1.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 0.5, 0.5); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 1.0, 0.8); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 2.0, 2.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 20.0, 20.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) dbeta(x, 3.0, 8.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) 0.9*dbeta(x, 3.0, 8.0) + 0.1*dbeta(x, 15, 1.5); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) 0.5*(x >= 0.5)*dbeta(x, 1.0, 2.0) + 0.5*(x < 0.5)*dbeta(x, 2.0, 1.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) 0.5*(x < 0.5)*(1.0 - 2.0*x) + 0.5*(x >= 0.5)*2.0*(x - 0.5); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) 0.4*dbeta(x, 30.0, 10.0) + 0.6*dbeta(x, 10.0, 30.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) 0.4*dbeta(x, 30.0, 3.0) + 0.6*dbeta(x, 3.0, 30.0) + 0.2; meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)
# h <- function(x) 2.0*dbeta(x, 30.0, 10.0) + 5.0*dbeta(x, 3, 3) + 1.7*dbeta(x, 10.0, 30.0); meanSliceWidth_int(h); meanSliceWidth_grid(xx, h(xx)); auc_w(h)


