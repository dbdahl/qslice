############## Independence Metropolis Hastings ################

IMH_sampler <- function(lf, x_0, pseudo_lpdf, pseudo_inv_cdf) {

  ## proposed draw
  x.dot <- pseudo_inv_cdf(runif(1, 0, 1))
  logr <- lf(x.dot) - lf(x_0) + pseudo_lpdf(x_0) - pseudo_lpdf(x.dot)
  u <- runif(1, 0, 1)
  
  if (log(u) < logr) {
    x_out <- x.dot
    accept <- 1
  } else {
    x_out <- x_0
    accept <- 0
  }

  list(x = x_out, accept = accept)
}


############## Random Walk ################

random_walk_sampler <- function(lf, support, x_0, c) {
  
  x_out <- x_0
  accept <- 0
  
  ## proposed draw
  x.dot <- rnorm(1, x_0, c)
  
  if (x.dot >= support[1] && x.dot <= support[2]) {
    
    logr <- lf(x.dot) - lf(x_0)
    u <- runif(1, 0, 1)
    
    if (log(u) < logr) {
      x_out <- x.dot
      accept <- 1
    }
  }
  
  list(x = x_out, accept = accept)
}
