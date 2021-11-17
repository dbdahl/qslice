#' Slice Sampler using the Stepping Out and Shrinkage Procedures
#'
#' This function implements a univariate slice sampler of Neal (2003) using the
#' "stepping out" procedure, followed by the "shrinkage" procedure.
#'
#' @param x The current state (as a numeric scalar).
#' @param target A function taking numeric scalar and returning a numeric
#'   scalar.
#' @param w A numeric scalar tuning the algorithm which gives the typical slice
#'   width. This is a main tuning parameter of the algorithm.
#' @param max The maximum number of times to step out. Setting \code{max} to
#'   zero avoids some evaluations of \code{target}, but may lead to relatively
#'   high autocorrelation if \code{w} is too small.  If \code{w} is too small,
#'   setting \code{max} to a large value (even \code{Inf}) should lead to low
#'   autocorrelation at the cost of more evaluations for \code{target}.
#' @param log If \code{FALSE}, the \code{target} function is the density
#'   function of the target distribution.  If \code{TRUE}, the function
#'   evaluates the log of the density.
#'
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log=TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' system.time({
#'     for ( i in seq.int(2,length(draws)) ) {
#'         out <- slice_sampler_stepping_out(draws[i-1], target=lf, w=0.7, max=Inf)
#'         draws[i] <- out$x
#'         nEvaluations <- nEvaluations + out$nEvaluations
#'     }
#' })
#' nEvaluations/length(draws)
#' plot(density(draws), xlim=c(0,1))
#' curve(exp(lf(x)), 0, 1, col="blue", add=TRUE)
#'
slice_sampler_stepping_out <- function(x, target, w, max=0, log=TRUE) {
  u <- function() runif(1, 0, 1)
  # Step 1
  nEvaluations <- 1
  fx <- target(x)
  y <- if ( isTRUE(log) ) {
    log(u() * exp(fx))
  } else {
    u() * fx
  }
  # Step 2 ("Stepping out" procedure)
  L <- x - u() * w
  R <- L + w
  if ( ! is.finite(max) ) {
    nEvaluations <- nEvaluations + 1
    while ( y < target(L) ) {
      nEvaluations <- nEvaluations + 1
      L <- L - w
    }
    nEvaluations <- nEvaluations + 1
    while ( y < target(R) ) {
      nEvaluations <- nEvaluations + 1
      R <- R + w
    }
  } else if ( max > 0 ) {
    J <- floor( u() * max )
    K <- max - 1 - J
    nEvaluations <- nEvaluations + 1
    while ( J > 0 && y < target(L) ) {
      nEvaluations <- nEvaluations + 1
      L <- L - w
      J <- J - 1
    }
    nEvaluations <- nEvaluations + 1
    while ( K > 0 && y < target(R) ) {
      nEvaluations <- nEvaluations + 1
      R <- R + w
      K <- K - 1
    }
  }
  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + u() * ( R - L )
    nEvaluations <- nEvaluations + 1
    fx1 <- target(x1)
    if ( y < fx1 ) return(list(x=x1, nEvaluations=nEvaluations))
    if ( x1 < x ) L <- x1 else R <- x1
  }
}

#' Transform Slice Sampler
#'
#' This function implements the transform slice sampler.  ...More details...
#'
#' @inherit slice_sampler_stepping_out
#' @param pseudo_log_pdf Not yet documented.
#' @param pseudo_log_cdf Not yet documented.
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log=TRUE)
#' pseudoLogPDF <- function(x) dbeta(x, shape1=1, shape2=1, log=TRUE)
#' pseudoInvCDF <- function(u) qbeta(u, shape1=1, shape2=1)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' system.time({
#'     for ( i in seq.int(2,length(draws)) ) {
#'         out <- slice_sampler_transform(draws[i-1], target=lf, pseudoLogPDF, pseudoInvCDF)
#'         draws[i] <- out$x
#'         nEvaluations <- nEvaluations + out$nEvaluations
#'     }
#' })
#' nEvaluations/length(draws)
#' plot(density(draws), xlim=c(0,1))
#' curve(exp(lf(x)), 0, 1, col="blue", add=TRUE)
#'
slice_sampler_transform <- function(x, target, pseudo_log_pdf, pseudo_inv_cdf, log=TRUE) {
  if ( ! isTRUE(log) ) stop("'log=FALSE' is not implemented.")
  u <- function() runif(1, 0, 1)
  # Step 1
  nEvaluations <- 1
  lfx <- target(x) - pseudo_log_pdf(x)
  y <- log(u() * exp(lfx))
  # Step 2 ("Shrinkage" procedure)
  L <- 0
  R <- 1
  repeat {
    u1 <- L + u() * ( R - L )
    x1 <- pseudo_inv_cdf(u1)
    nEvaluations <- nEvaluations + 1
    lfx1 <- target(x1) - pseudo_log_pdf(x1)
    if ( y < lfx1 ) return(list(x=x1, nEvaluations=nEvaluations))
    if ( x1 < x ) L <- u1 else R <- u1
  }
}

#' Latent Slice Sampler
#'
#' This function implements the latent slice sampler of Li and Walker (2020).
#' ...More details...
#'
#' @inherit slice_sampler_stepping_out
#' @param s A random variable that determines how far the algorithm samples from
#'   on each side
#' @param rate The rate parameter for a truncated exponential.
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log=TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' system.time({
#'     for ( i in seq.int(2,length(draws)) ) {
#'         out <- slice_sampler_latent(draws[i-1], 0.5, target=lf, rate=1)
#'         draws[i] <- out$x
#'         nEvaluations <- nEvaluations + out$nEvaluations
#'     }
#' })
#' nEvaluations/length(draws)
#' plot(density(draws), xlim=c(0,1))
#' curve(exp(lf(x)), 0, 1, col="blue", add=TRUE)
#'
slice_sampler_latent <- function(x, s, target, rate, log=TRUE) {
  if ( ! isTRUE(log) ) stop("'log=FALSE' is not implemented.")
  nEvaluations <- 1
  lfx <- target(x)
  ly <- log(runif(1)) + lfx
  half_s <- s/2
  l <- runif(1, x - half_s, x + half_s)
  # Eq. 7... a truncated exponential using the inverse CDF method.
  s <- -log(runif(1))/rate + 2*abs(l-x)
  half_s <- s/2
  L <- l - half_s
  R <- l + half_s
  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- L + runif(1) * ( R - L )
    nEvaluations <- nEvaluations + 1
    lfx1 <- target(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, s=s, nEvaluations=nEvaluations)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}

#' Univariate Elliptical Slice Sampler
#'
#' This function implements the elliptical slice sampler of Murray, Adams, MacKay (2010)
#' ...More details...
#'
#' @inherit slice_sampler_stepping_out
#' @param mu A numeric scalar tuning the algorithm which gives the theta value that will be used to sample a random value from the ellipse.
#' @param sigma A numeric scalar tuning the algorithm which gives the theta value that will be used to sample a random value from the ellipse.
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log=TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' system.time({
#'     for ( i in seq.int(2,length(draws)) ) {
#'         out <- slice_sampler_elliptical(draws[i-1], target=lf, mu=0.5, sigma=1)
#'         draws[i] <- out$x
#'         nEvaluations <- nEvaluations + out$nEvaluations
#'     }
#' })
#' nEvaluations/length(draws)
#' plot(density(draws), xlim=c(0,1))
#' curve(exp(lf(x)), 0, 1, col="blue", add=TRUE)
#'
slice_sampler_elliptical <- function(x, target, mu = 2, sigma = 5) {
  nEvaluations <- 1
  nu <- rnorm(1,mu,sigma)
  u <- runif(1,0,1)
  log_y <- log(target(x)) + log(u)
  theta <- runif(1,0,2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  # Sample x1
  x1 <- (x - mu)*cos(theta) + (nu - mu)*sin(theta) + mu
  # Checking to see if x1 is in the distribution. If not then shrinkage procedure for theta
  while (log(target(x1)) < log_y) {
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
  results
}

#' General Elliptical Slice Sampler
#'
#' General Elliptical Slice Sampler of Nishihara (2014)
#'
#' @inheritParams slice_sampler_stepping_out
#' @inheritParams slice_sampler_elliptical
#' @param degf Number of degrees of freedom
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log=TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' system.time({
#'     for ( i in seq.int(2,length(draws)) ) {
#'         out <- slice_sampler_generalized_elliptical(draws[i-1], target=lf, mu=0.5, sigma=1, degf=5)
#'         draws[i] <- out$x
#'         nEvaluations <- nEvaluations + out$nEvaluations
#'     }
#' })
#' nEvaluations/length(draws)
#' plot(density(draws), xlim=c(0,1))
#' curve(exp(lf(x)), 0, 1, col="blue", add=TRUE)
#'
slice_sampler_generalized_ellipse <- function(x, target, mu, sigma, degf) {
  ## here, f is the target density
  nEvaluations <- 1
  a <- (degf + 1.0) / 2.0
  b <- 0.5*(degf + ((x - mu)/sigma)^2)
  s <- 1.0 / rgamma(1, shape=a, rate=b) # rate of gamma <=> shape of inv-gamma
  lff <- function(xx) target(xx) - ( dt((xx-mu)/sigma, df=degf, log=TRUE) - log(sigma) )
  slice_sampler_elliptical(x=x, target=lff, mu=mu, sigma=sqrt(s)*sigma)
}

# #' Slice Sampler using the Stepping Out and Shrinkage Procedures (using Rust)
# #'
# #' @inheritParams slice_sampler_stepping_out
# #' @inherit slice_sampler_stepping_out return examples
# #'
# slice_sampler <- function(x, target, w, max=0, log=FALSE) {
#   .Call(.slice_sampler, x, target, w, max, log)
# }
