#' A function that takes a positive function and returns a vector
#' @param lx any function that is always positive
#' @param mu the mean of the normal distribution you will be sampling from
#' @param sigma the standard deviation of the normal distribution you will be sampling from
#' @param iters the number of iterations you want to sample
#' @param start_x the starting value of x
#' @return returns a vector of x values in your sample space
#' @import stats
#' @export

elliptical_slice_sample <- function(lx, mu = .2, sigma = 100, iters = 1000, start_x = 0){
  x_0 <- vector(mode = 'numeric', length = iters)
  y <- vector(mode = 'numeric', length = iters)
  nEvaluations <- 1
  # Staring the loop for number of iterations
  for(i in 1:iters) {
    if (i == 1) {
      x_0[i] <- start_x
    } else {
      x_0[i] <- x_1
    }
    v <- rnorm(1,mu,sigma)
    u <- runif(1,0,1)
    log_y <- log(lx(x_0[i])) + log(u)
    y[i] <- exp(log_y)
    theta <- runif(1,0,2*pi)
    theta_min <- theta - 2*pi
    theta_max <- theta
    x_1 <- (x_0[i] - mu)*cos(theta) + (v - mu)*sin(theta) + mu
    # Shrinkage procedure
    while (log(lx(x_1)) < log_y) {
      if (theta < 0){
        theta_min <- theta
      } else {
        theta_max <- theta
      }
      theta <- runif(1,theta_min, theta_max)
      nEvaluations <- nEvaluations + 1
      x_1 <- (x_0[i] - mu)*cos(theta) + (v - mu)*sin(theta) + mu
    }
    x_0[i] <- x_1
  }
  sample_space <- data.frame("x" = x_0, "y" = y)
  return(list(sample_space, nEvaluations))
}

