#' Slice sampler using the Stepping Out and Shrinkage Procedures
#'
#' Single update for the univariate slice sampler of Neal (2003) using the
#' "stepping out" procedure, followed by the "shrinkage" procedure.
#'
#' @param x The current state (as a numeric scalar).
#' @param log_target A function taking numeric scalar that evaluates the
#' (potentially unnormalized) log-target density, returning a numeric scalar.
#' @param w A numeric scalar tuning the algorithm which gives the typical slice
#'   width. This is a main tuning parameter of the algorithm.
#' @param max The maximum number of times to step out. Setting \code{max} to
#'   zero avoids some evaluations of \code{log_target}, but may lead to relatively
#'   high autocorrelation if \code{w} is too small.  If \code{w} is too small,
#'   setting \code{max} to a large value (even \code{Inf}) should lead to low
#'   autocorrelation at the cost of more evaluations for \code{log_target}.
#'
#' @return A list with two elements:
#'
#' \code{x} is the new state.
#'
#' \code{nEvaluations} is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @references
#' Neal, R. M. (2003), "Slice sampling," *The Annals of Statistics*, 31, 705-767. \doi{https://doi.org/10.1214/aos/1056562461}
#'
#' @importFrom stats runif
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(10) + 0.5 # set to numeric(1e3) for more complete illustration
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_stepping_out(draws[i - 1], log_target = lf, w = 0.7, max = Inf)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (length(draws) - 1)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_stepping_out <- function(x, log_target, w, max = Inf) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x)
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

#' Quantile Slice Sampler
#'
#' Single update using a quantile slice sampler of Heiner et al. (2024+).
#'
#' @inherit slice_stepping_out
#' @param pseudo List containing two functions specifying the pseudo-target distribution:
#'
#' \code{ld} evaluates the log density for a scalar input, and
#'
#' \code{q} evaluates the quantile (inverse-CDF) function for an input in (0,1).
#'
#' @return A list containing three elements:
#'
#' \code{x} is the new state.
#'
#' \code{u} is the value of the CDF of the psuedo-target associated with the
#' returned value (also referred to as psi).
#'
#' \code{nEvaluations} is the number of evaluations of the
#'   target function used to obtain the new state.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*.
#'
#' @importFrom stats runif
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' pseu <- list(ld = function(x) dbeta(x, shape1 = 1, shape2 = 1, log = TRUE),
#'              q = function(u) qbeta(u, shape1 = 1, shape2 = 1))
#' draws <- numeric(10) # set to numeric(1e3) for more complete illustration
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_quantile(draws[i - 1], log_target = lf, pseudo = pseu)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (length(draws) - 1)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_quantile <- function(x, log_target, pseudo) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x) - pseudo$ld(x)
  }
  # Step 1
  y <- log(runif(1)) + f(x)
  # Step 2 ("Shrinkage" procedure)
  L <- 0
  R <- 1
  repeat {
    u1 <- runif(1, L, R)
    x1 <- pseudo$q(u1)
    if (y < f(x1)) {
      return(list(x = x1, u = u1, nEvaluations = nEvaluations))
    }
    if (x1 < x) L <- u1 else R <- u1
  }
}

#' Latent Slice Sampler
#'
#' Single update using the latent slice sampler of Li and Walker (2023).
#'
#' @inherit slice_stepping_out
#' @param s A random variable that determines the length of the initial shrinking interval.
#' @param rate The rate parameter for the distribution of \code{s}.
#'
#' @return A list containing three elements:
#'
#' \code{x} is the new state of the target variable.
#'
#' \code{s} is the new state of the latent scale variable.
#'
#' \code{nEvaluations} is the number of evaluations of the
#'   target function used to obtain the new state.
#'
#' @references
#' Li, Y. and Walker, S. G. (2023), "A latent slice sampling algorithm," *Computational Statistics and Data Analysis*, 179, 107652. \doi{https://doi.org/10.1016/j.csda.2022.107652}
#'
#' @importFrom stats runif
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(10) # set to numeric(1e3) for more complete illustration
#' nEvaluations <- 0L
#' s <- 0.5
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_latent(draws[i - 1], s, log_target = lf, rate = 0.3)
#'   draws[i] <- out$x
#'   s <- out$s
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (length(draws) - 1)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_latent <- function(x, s, log_target, rate) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x)
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
#' Algorithm 1 of Nishihara et al. (2014) of the
#' elliptical slice sampler of Murray et al. (2010).
#'
#' @inherit slice_stepping_out
#' @param mu A numeric scalar with the mean of the supporting normal distribution.
#' @param sigma A numeric scalar with the standard deviation of the supporting normal distribution.
#'
#' @references
#' Murray, I., Adams, R., and MacKay, D., (2010), "Elliptical Slice Sampling," in *Proceedings of the Thirteenth International Conference on Artificial Intelligence and Statistics*, JMLR Workshop and Conference Proceedings. \url{https://proceedings.mlr.press/v9/murray10a}
#'
#' Nishihara, R., Murray, I., and Adams, R. P. (2014), "Parallel MCMC with Generalized Elliptical Slice Sampling," *Journal of Machine Learning Research*, 15, 2087-2112. \url{https://jmlr.org/papers/v15/nishihara14a.html}
#'
#' @importFrom stats runif rnorm
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(10) # set to numeric(1e3) for more complete illustration
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_elliptical(draws[i - 1], log_target = lf, mu = 0.5, sigma = 1)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (length(draws) - 1)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_elliptical <- function(x, log_target, mu, sigma) {
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x)
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

#' Generalized Elliptical Slice Sampler (univariate)
#'
#' Single update using the generalized elliptical slice sampler of Nishihara et al. (2014).
#'
#' @inheritParams slice_stepping_out
#' @inheritParams slice_elliptical
#' @param df Degrees of freedom of Student t pseudo-target.
#'
#' @return A list contains two elements:
#'
#' \code{x} is the new state.
#'
#' \code{nEvaluations} is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @references
#' Nishihara, R., Murray, I., and Adams, R. P. (2014), "Parallel MCMC with Generalized Elliptical Slice Sampling," *Journal of Machine Learning Research*, 15, 2087-2112. \url{https://jmlr.org/papers/v15/nishihara14a.html}
#'
#' @importFrom stats dt rgamma
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(10) # set to numeric(1e3) for more complete illustration
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_genelliptical(draws[i - 1], log_target = lf,
#'                                       mu = 0.5, sigma = 1, df = 5)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (length(draws) - 1)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
slice_genelliptical <- function(x, log_target, mu, sigma, df) {
  a <- (df + 1.0) / 2.0
  b <- 0.5 * (df + ((x - mu) / sigma)^2)
  s <- 1.0 / rgamma(1, shape = a, rate = b) # rate of gamma <=> shape of inv-gamma
  lff <- function(xx) log_target(xx) - (dt((xx - mu) / sigma, df = df, log = TRUE) - log(sigma))
  slice_elliptical(x = x, log_target = lff, mu = mu, sigma = sqrt(s) * sigma)
}
