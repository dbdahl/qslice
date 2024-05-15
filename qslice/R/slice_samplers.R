#' Slice Sampler using the Stepping Out and Shrinkage Procedures
#'
#' Univariate slice sampler of Neal (2003) using the
#' "stepping out" procedure, followed by the "shrinkage" procedure.
#'
#' @param x The current state (as a numeric scalar).
#' @param target A function taking numeric scalar that evaluates the log-target
#'   density, returning a numeric scalar.
#' @param w A numeric scalar tuning the algorithm which gives the typical slice
#'   width. This is a main tuning parameter of the algorithm.
#' @param max The maximum number of times to step out. Setting \code{max} to
#'   zero avoids some evaluations of \code{target}, but may lead to relatively
#'   high autocorrelation if \code{w} is too small.  If \code{w} is too small,
#'   setting \code{max} to a large value (even \code{Inf}) should lead to low
#'   autocorrelation at the cost of more evaluations for \code{target}.
#'
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_stepping_out(draws[i - 1], target = lf, w = 0.7, max = Inf)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / length(draws)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_stepping_out <- function(x, target, w, max = Inf) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
  }
  # Step 1
  fx <- f(x)
  y <- log(runif(1)) + fx
  # Step 2 ("Stepping out" procedure)
  L <- x - runif(1) * w
  R <- L + w
  if (!is.finite(max)) {
    while (y < f(L)) L <- L - w
    while (y < f(R)) R <- R + w
  } else if (max > 0) {
    J <- floor(runif(1) * max)
    K <- max - 1 - J
    while (J > 0 && y < f(L)) {
      L <- L - w
      J <- J - 1
    }
    while (K > 0 && y < f(R)) {
      R <- R + w
      K <- K - 1
    }
  }
  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + runif(1) * (R - L)
    if (y < f(x1)) {
      return(list(x = x1, nEvaluations = nEvaluations))
    }
    if (x1 < x) L <- x1 else R <- x1
  }
}

#' Transform Slice Sampler
#'
#' Quantile slice sampler.
#'
#' @inherit slice_stepping_out
#' @param pseudo_log_pdf Not yet documented.
#' @param pseudo_inv_cdf Not yet documented.
#'
#' @return A list containing three elements: "x" is the new state, "u" is the
#'   value of the CDF of the psuedo-target associated with the returned value,
#'   inverse CDF method, and "nEvaluations is the number of evaluations of the
#'   target function used to obtain the new state.
#'
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' pseudoLogPDF <- function(x) dbeta(x, shape1 = 1, shape2 = 1, log = TRUE)
#' pseudoInvCDF <- function(u) qbeta(u, shape1 = 1, shape2 = 1)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_transform(draws[i - 1], target = lf, pseudoLogPDF, pseudoInvCDF)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / length(draws)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_transform <- function(x, target, pseudo_log_pdf, pseudo_inv_cdf) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x) - pseudo_log_pdf(x)
  }
  # Step 1
  y <- log(runif(1)) + f(x)
  # Step 2 ("Shrinkage" procedure)
  L <- 0
  R <- 1
  repeat {
    u1 <- runif(1, L, R)
    x1 <- pseudo_inv_cdf(u1)
    if (y < f(x1)) {
      return(list(x = x1, u = u1, nEvaluations = nEvaluations))
    }
    if (x1 < x) L <- u1 else R <- u1
  }
}

#' Latent Slice Sampler
#'
#' Latent slice sampler of Li and Walker (2020).
#'
#' @inherit slice_stepping_out
#' @param s A random variable that determines how far the algorithm samples from
#'   on each side
#' @param rate The rate parameter for a truncated exponential.
#'
#' @return A list contains three elements: "x" is the new state of the target variable,
#'   "s" is the new state of the latent scale variable, and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' s <- 0.5
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_latent(draws[i - 1], s, target = lf, rate = 0.3)
#'   draws[i] <- out$x
#'   s <- out$s
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / length(draws)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_latent <- function(x, s, target, rate) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
  }
  # Step 1
  y <- log(runif(1)) + f(x)
  half_s <- s / 2
  l <- runif(1, x - half_s, x + half_s)
  # Eq. 7... a truncated exponential using the inverse CDF method.
  s <- -log(runif(1)) / rate + 2 * abs(l - x)
  half_s <- s / 2
  L <- l - half_s
  R <- l + half_s
  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- L + runif(1) * (R - L)
    if (y < f(x1)) {
      return(list(x = x1, s = s, nEvaluations = nEvaluations))
    }
    if (x1 < x) L <- x1 else R <- x1
  }
}

#' Univariate Elliptical Slice Sampler
#'
#' Algorithm 1 of Nishihara et al (2014) of the
#' elliptical slice sampler of Murray, Adams, MacKay (2010).
#'
#' @inherit slice_stepping_out
#' @param mu A numeric scalar with the mean of the normal to be sampled
#' @param sigma A numeric scalar with the standard deviation of the normal to be sampled.
#'
#' @importFrom stats runif rnorm
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_elliptical(draws[i - 1], target = lf, mu = 0.5, sigma = 1)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / length(draws)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_elliptical <- function(x, target, mu, sigma) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
  }
  # Step 1
  y <- log(runif(1)) + f(x)
  nu <- rnorm(1, mu, sigma)
  theta <- runif(1, 0, 2 * pi)
  theta_min <- theta - 2 * pi
  theta_max <- theta
  repeat {
    x1 <- (x - mu) * cos(theta) + (nu - mu) * sin(theta) + mu
    if (y < f(x1)) {
      return(list(x = x1, nEvaluations = nEvaluations))
    }
    if (theta < 0) {
      theta_min <- theta
    } else {
      theta_max <- theta
    }
    theta <- runif(1, theta_min, theta_max)
  }
}

#' General Elliptical Slice Sampler (Univariate)
#'
#' General Elliptical Slice Sampler of Nishihara (2014)
#'
#' @inheritParams slice_stepping_out
#' @inheritParams slice_elliptical
#' @param df Degrees of freedom of Student t pseudo-target.
#'
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @importFrom stats dt rgamma
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(1000)
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_generalized_elliptical(draws[i - 1], target = lf,
#'                                       mu = 0.5, sigma = 1, df = 5)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / length(draws)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_generalized_elliptical <- function(x, target, mu, sigma, df) {
  a <- (df + 1.0) / 2.0
  b <- 0.5 * (df + ((x - mu) / sigma)^2)
  s <- 1.0 / rgamma(1, shape = a, rate = b) # rate of gamma <=> shape of inv-gamma
  lff <- function(xx) target(xx) - (dt((xx - mu) / sigma, df = df, log = TRUE) - log(sigma))
  slice_elliptical(x = x, target = lff, mu = mu, sigma = sqrt(s) * sigma)
}

