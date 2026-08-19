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
  stopifnot(
    length(x) == 1L,
    is.finite(x),
    is.function(log_target),
    length(w) == 1L,
    is.finite(w),
    w > 0.0,
    is.numeric(max),
    length(max) == 1L,
    !is.na(max),
    (is.infinite(max) && max > 0) ||
      (is.finite(max) && max >= 0 && max == floor(max))
  )

  nEvaluations <- 0L

  lf <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_target,
      x = z,
      sampler = "slice_stepping_out",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }

  # Step 1
  lfx <- lf(x, stage = "evaluate at current state", allow_minus_inf = FALSE)
  y <- log(runif(1, min = 0.0, max = 1.0)) + lfx
  # Step 2 ("Stepping out" procedure)
  L <- x - runif(1, min = 0.0, max = w)
  R <- L + w
  if (!is.finite(max)) {
    while (y < lf(L, stage = "left stepping-out proposal")) {
      L <- L - w
    }
    while (y < lf(R, stage = "right stepping-out proposal")) {
      R <- R + w
    }
  } else if (max > 0) {
    J <- floor(runif(1, min = 0.0, max = max))
    K <- max - 1 - J
    while (J > 0 && y < lf(L, stage = "left stepping-out proposal")) {
      L <- L - w
      J <- J - 1
    }
    while (K > 0 && y < lf(R, stage = "right stepping-out proposal")) {
      R <- R + w
      K <- K - 1
    }
  }
  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- runif(1, min = L, max = R)
    if (y < lf(x1, stage = "shrinkage proposal")) {
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
#' \code{u} is the value of the CDF of the pseudo-target associated with the
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
#' draws <- numeric(10) + 0.5 # set to numeric(1e3) + 0.5 for more complete illustration
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
  stopifnot(
    length(x) == 1L,
    is.finite(x),
    is.function(log_target),
    is.list(pseudo),
    is.function(pseudo$ld),
    is.function(pseudo$q)
  )

  nEvaluations <- 0L

  lf <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_target,
      x = z,
      sampler = "slice_quantile",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }

  ld_pseudo <- function(z, stage) {
    .eval_pseudo_logdens(
      logdens = pseudo$ld,
      x = z,
      sampler = "slice_quantile",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = FALSE
    )
  }

  lh <- function(x, stage, allow_minus_inf = TRUE) {
    lf(x, stage = stage, allow_minus_inf = allow_minus_inf) -
      ld_pseudo(x, stage = stage)
  }

  # Step 1
  y <- log(runif(1, min = 0.0, max = 1.0)) +
    lh(x, stage = "evaluate at current state", allow_minus_inf = FALSE)
  # Step 2 ("Shrinkage" procedure)
  L <- 0.0
  R <- 1.0
  proposal_attempts <- 0L
  repeat {
    proposal_attempts <- proposal_attempts + 1L
    u1 <- runif(1, min = L, max = R)
    x1 <- .eval_pseudo_quantile(
      q = pseudo$q,
      u = u1,
      sampler = "slice_quantile",
      stage = "shrinkage proposal",
      attempt = proposal_attempts,
      expected_length = 1L
    )
    if (y < lh(x1, stage = "shrinkage")) {
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
#' draws <- numeric(10) + 0.5 # set to numeric(1e3) + 0.5 for more complete illustration
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
  stopifnot(
    length(x) == 1L,
    length(s) == 1L,
    length(rate) == 1L,
    is.finite(x),
    is.function(log_target),
    is.finite(s),
    is.finite(rate),
    s > 0.0,
    rate > 0.0
  )

  nEvaluations <- 0L
  lf <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_target,
      x = z,
      sampler = "slice_latent",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }
  # Step 1
  y <- log(runif(1, min = 0.0, max = 1.0)) +
    lf(x, stage = "evaluate at current state", allow_minus_inf = FALSE)
  half_s <- s / 2
  l <- runif(1, min = x - half_s, max = x + half_s)
  # Eq. 7... a truncated exponential using the inverse CDF method.
  s <- -log(runif(1, min = 0.0, max = 1.0)) / rate + 2.0 * abs(l - x)
  half_s <- s / 2.0
  L <- l - half_s
  R <- l + half_s
  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- runif(1, min = L, max = R)
    if (y < lf(x1, stage = "shrinkage step")) {
      return(list(x = x1, s = s, nEvaluations = nEvaluations))
    }
    if (x1 < x) L <- x1 else R <- x1
  }
}

.slice_elliptical_impl <- function(x, log_lik, mu, sigma, sampler) {
  nEvaluations <- 0L
  llik <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_lik,
      x = z,
      sampler = sampler,
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }

  # Step 1
  y <- log(runif(1, min = 0.0, max = 1.0)) +
    llik(x, stage = "evaluate at current state", allow_minus_inf = FALSE)
  nu <- rnorm(1, mean = mu, sd = sigma)
  theta <- runif(1, 0, 2 * pi)
  theta_min <- theta - 2 * pi
  theta_max <- theta
  repeat {
    x1 <- (x - mu) * cos(theta) + (nu - mu) * sin(theta) + mu
    if (y < llik(x1, stage = "shrinkage step")) {
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

#' Univariate Elliptical Slice Sampler
#'
#' Algorithm 1 of Nishihara et al. (2014) of the
#' elliptical slice sampler of Murray et al. (2010).
#'
#' @param x The current state (as a numeric scalar).
#' @param log_lik A function taking numeric scalar that evaluates the natural logarithm of the
#' "likelihood" L(x) part of the unnormalized target density L(x) * Normal(x; mu, sigma). Returns a numeric scalar.
#' @param mu A numeric scalar with the mean of the supporting normal "prior" distribution.
#' @param sigma A numeric scalar with the standard deviation of the supporting normal "prior" distribution.
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
#' log_lik <- function(x) dbeta(x, 8, 2, log = TRUE)
#' lf <- function(x) log_lik(x) + dnorm(x, mean = mu, sd = sig, log = TRUE)
#' mu <- 0.2
#' sig <- 0.3
#' draws <- numeric(10) + 0.5 # set to numeric(1e3) + 0.5 for more complete illustration
#' nEvaluations <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- slice_elliptical(draws[i - 1], log_lik = log_lik, mu = mu, sigma = sig)
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (length(draws) - 1)
#' plot(density(draws), xlim = c(-0.2, 1.2))
#' curve(exp(lf(x))*4.5, -0.2, 1.2, col = "blue", add = TRUE) # 4.5 approximates the normalizing constant
#'
slice_elliptical <- function(x, log_lik, mu, sigma) {
  stopifnot(
    length(x) == 1L,
    length(mu) == 1L,
    length(sigma) == 1L,
    is.finite(x),
    is.function(log_lik),
    is.finite(mu),
    is.finite(sigma),
    sigma > 0.0
  )

  .slice_elliptical_impl(
    x = x,
    log_lik = log_lik,
    mu = mu,
    sigma = sigma,
    sampler = "slice_elliptical"
  )
}

#' Generalized Elliptical Slice Sampler (univariate)
#'
#' Single update using the generalized elliptical slice sampler of Nishihara et al. (2014).
#'
#' @inheritParams slice_stepping_out
#' @param mu A numeric scalar with the location parameter of the Student t pseudo-target.
#' @param sigma A positive numeric scalar with the scale parameter of the Student t pseudo-target.
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
#' draws <- numeric(10) + 0.5 # set to numeric(1e3) + 0.5 for more complete illustration
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
  stopifnot(
    length(x) == 1L,
    length(mu) == 1L,
    length(sigma) == 1L,
    length(df) == 1L,
    is.finite(x),
    is.function(log_target),
    is.finite(mu),
    is.finite(sigma),
    is.finite(df),
    sigma > 0.0,
    df > 0.0
  )

  a <- (df + 1.0) / 2.0
  b <- 0.5 * (df + ((x - mu) / sigma)^2)
  s <- 1.0 / rgamma(1, shape = a, rate = b) # rate of gamma <=> shape of inv-gamma
  lff <- function(xx) {
    log_target(xx) - (dt((xx - mu) / sigma, df = df, log = TRUE)) #  - log(sigma) is const. and not necessary
  }
  .slice_elliptical_impl(
    x = x,
    log_lik = lff,
    mu = mu,
    sigma = sqrt(s) * sigma,
    sampler = "slice_genelliptical"
  )
}
