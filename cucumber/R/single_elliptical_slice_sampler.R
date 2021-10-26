#' Slice Sampler using the Stepping Out and Shrinkage Procedures
#'
#' This function implements a univariate elliptical slice sampler of Nishihara (2014)
#'
#' @param x The current state (as a numeric scalar).
#' @param fx A function taking numeric scalar and returning a numeric
#'   scalar.
#' @param mu A numeric scalar tuning the algorithm which gives the theta value that will be used to sample a random value from the ellipse
#' @param sigma A numeric scalar tuning the algorithm which gives the theta value that will be used to sample a random value from the ellipse
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @import stats
#' @export
single_elliptical_slice_sample <- function(x = 0, mu = 2, sigma = 5, fx) {
  nEvaluations <- 1
  v <- rnorm(1,mu,sigma)
  u <- runif(1,0,1)
  log_y <- log(fx(x)) + log(u)
  theta <- runif(1,0,2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  # Sample x1
  x1 <- (x - mu)*cos(theta) + (v - mu)*sin(theta) + mu
  # Checking to see if x1 is in the distribution. If not then shrinkage procedure
  while (log(fx(x1)) < log_y) {
    if (theta < 0){
      theta_min <- theta
    } else {
      theta_max <- theta
    }
    theta <- runif(1,theta_min, theta_max)
    nEvaluations <- nEvaluations + 1
    x1 <- (x - mu)*cos(theta) + (v - mu)*sin(theta) + mu
  }
  results <- list(x = x1, nEvaluations = nEvaluations)
  return(results)
}
