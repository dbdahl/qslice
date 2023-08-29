skew_hist <- function(x, y, u = NULL, nbins = 30) {
  
  if (is.null(u)) {
    
    if (is.function(y)) {
      y <- y(x)
    } 
    
    stopifnot (length(y) == length(x))
    stopifnot(all(x > 0.0) && all(y >= 0.0) && all(x < 1.0))
    
  } else { ## a less noisy estimate would be to calculate sample skewness
    
    bins <- seq(0.0, 1.0, len = nbins + 1)
    y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins=nbins)
    x <- (bins[-(nbins+1)] + bins[-1]) / 2.0
    
  }
  
  yn <- y / sum(y) # distribution rather than modal value of 1
  
  mu <- sum(x * yn)
  mu_2 <- sum(x^2 * yn)
  mu_3 <- sum(x^3 * yn)
  
  mu2 <- mu^2
  mu3 <- mu^3
  
  sig2 <- mu_2 - mu2
  sig <- sqrt(sig2)
  
  ( mu_3 - 3.0*mu*sig2 - mu3 ) / sig^3
  
}

# n = 30
# x <- seq(0.01, 0.99, len = n)
# 
# y <- dbeta(x, 4, 2)
# y <- dbeta(x, 1, 2)
# y <- dbeta(x, 2, 1)
# y <- dbeta(x, 0.5, 0.5)
# y <- dbeta(x, 0.2, 0.8)
# y <- dbeta(x, 1, 1)
# 
# plot(x, y)
# 
# skew_hist(x, y)
# 
# 
# library("moments")
# skewness(rbeta(10e3, 0.2, 0.8))

