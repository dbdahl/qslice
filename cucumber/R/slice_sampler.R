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
slice_sampler_stepping_out <- function(x, target, w, max=0, log=FALSE) {
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
slice_sampler_latent <- function(x, s, lf, rate) {
  nEvaluations <- 1
  lfx <- lf(x)
  ly <- log(stunif()) + lfx
  half_s <- s/2
  l <- runif(1, x - half_s, x + half_s)
  # Eq. 7... a truncated exponential using the inverse CDF method.
  s <- -log(runif(1))/rate + 2*abs(l-x)
  half_s <- s/2
  L <- l - half_s
  R <- l + half_s
  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    nEvaluations <- nEvaluations + 1
    lfx1 <- lf(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, s=s, nEvaluations=nEvaluations)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}
