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
#'   function of target distribution.  If \code{TRUE}, the function is log of
#'   the density.
#'
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @export
#' @examples
#' f <- function(x) dbeta(x, 3, 4, log=TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' system.time({
#' for ( i in seq_along(draws)[-1] ) {
#'     out <- slice_sampler(draws[i-1], f, w=0.7, max=Inf, log=TRUE)
#'     draws[i] <- out$x
#'     nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' })
#' nEvaluations/length(draws)
#' plot(density(draws), xlim=c(0,1))
#' curve(exp(f(x)), 0, 1, col="blue", add=TRUE)
#'
slice_sampler <- function(x, target, w, max=0, log=FALSE) {
  .Call(.slice_sampler, x, target, w, max, log)
}

#' @export
slice_sampler_ <- function(x, target, w, max=0, log=FALSE) {
  u <- function() runif(1, 0, 1)
  # Step 1
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
    while ( y < target(L) ) {
      L <- L - w
    }
    while ( y < target(R) ) {
      R <- R + w
    }
  } else if ( max > 0 ) {
    J <- floor( u() * max )
    K <- max - 1 - J
    while ( J > 0 && y < target(L) ) {
      L <- L - w
      J <- J - 1
    }
    while ( K > 0 && y < target(R) ) {
      R <- R + w
      K <- K - 1
    }
  }
  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + u() * ( R - L )
    fx1 <- target(x1)
    if ( y < fx1 ) return(x1)
    if ( x1 < x ) L <- x1 else R <- x1
  }
}

#' Transform Slice Sampler
#'
#' @export
slice_sampler_transform <- function(x, log_density, pseudo_log_pdf, pseudo_inv_cdf) {
  u <- function() runif(1, 0, 1)
  # Step 1
  lfx <- log_density(x) - pseudo_log_pdf(x)
  y <- log(u() * exp(lfx))
  # Step 2 ("Shrinkage" procedure)
  L <- 0
  R <- 1
  repeat {
    u1 <- L + stunif() * ( R - L )
    x1 <- pseudo_inv_cdf(u1)
    lfx1 <- log_density(x1) - pseudo_log_pdf(x1)
    if ( y < lfx1 ) return(x1)
    if ( x1 < x ) L <- u1 else R <- u1
  }
}

#' Algorithm of Li and Walker (2020)
#'
#' @export
slice_sampler_latent <- function(x, f, lambda) {

  # Step 1
  if ( is.list(x) ) {
    fx <- x$fx
    s <- x$s
    x <- x$x
  } else {
    stop("Supplied x must be a list with previous value, previous latent s, and density eval at previous value.")
  }

  y <- stunif() * fx

  # Step 2 sample latent variates
  ell <- runif(1, min = x - 0.5*s, max = x + 0.5*s)

  s_trunc <- 2.0*abs(ell - x)
  qtile_trunc <- 1.0 - exp(-lambda*s_trunc)
  s1 <- -log(1.0 - runif(1, min=qtile_trunc, max=1.0)) / lambda # draw from truncated exponential

  L <- ell - 0.5*s1
  R <- ell + 0.5*s1

  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    fx1 <- f(x1)
    if ( y < fx1 ) {
      result <- list(x=x1, s=s1, fx=fx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}
