.slice_elliptical_mv_impl <- function(
  x,
  log_lik,
  mu,
  Sig,
  is_chol = FALSE,
  sampler
) {
  K <- length(x)
  stopifnot(
    length(mu) == K,
    dim(Sig) == c(K, K),
    all(is.finite(x)),
    is.function(log_lik),
    all(is.finite(mu)),
    all(is.finite(Sig)),
    is.logical(is_chol)
  )

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

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
  llik_x0 <- llik(
    x,
    stage = "evaluate at current state",
    allow_minus_inf = FALSE
  )
  y <- log(runif(1, min = 0.0, max = 1.0)) + llik_x0

  nu <- drop(SigL %*% rnorm(K, mean = 0.0, sd = 1.0) + mu)

  twopi <- 2.0 * pi
  theta <- runif(1, min = 0, max = twopi)
  theta_min <- theta - twopi
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
    theta <- runif(1, min = theta_min, max = theta_max)
  }
}

#' Multivariate Elliptical Slice Sampler
#'
#' Algorithm 1 of Nishihara et al. (2014) of the
#' elliptical slice sampler of Murray et al. (2010).
#'
#' @param x The current state (as a numeric vector).
#' @param log_lik A function taking numeric vector that evaluates the natural logarithm of the
#' "likelihood" L(x) part of the unnormalized target density L(x) * Normal(x; mu, Sig). Returns a numeric scalar.
#' @param mu Numeric vector with the mean of the supporting normal "prior" distribution.
#' @param Sig Positive definite covariance matrix of the supporting normal "prior" distribution. Alternatively, a
#' lower-triangular matrix with the Cholesky factor of the covariance matrix
#' (for faster computation).
#' @param is_chol Logical, is the supplied \code{Sig} in Cholesky (lower triangular) format? Default is false.
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
#' mu <- c(0.5, 0.5)
#' Sig <- matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2)
#' SigL <- t(chol(Sig))
#' log_lik <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' lf <- function(x) log_lik(x) - 0.5 * drop(crossprod(forwardsolve(SigL, x - mu)))
#' n_iter <- 10 # set to 1e3 for more complete illustration
#' draws <- matrix(0.3, nrow = n_iter, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, n_iter)) {
#'   out <- slice_elliptical_mv(draws[i - 1,], log_lik = log_lik,
#'               mu = mu, Sig = SigL, is_chol = TRUE)
#'   draws[i,] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (n_iter - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#'
slice_elliptical_mv <- function(
  x,
  log_lik,
  mu,
  Sig,
  is_chol = FALSE
) {
  .slice_elliptical_mv_impl(
    x = x,
    log_lik = log_lik,
    mu = mu,
    Sig = Sig,
    is_chol = is_chol,
    sampler = "slice_elliptical_mv"
  )
}

#' Generalized Elliptical Slice Sampler (Multivariate)
#'
#' Generalized Elliptical Slice Sampler, Algorithm 2 of Nishihara et al. (2014)
#'
#' @param x The current state (as a numeric vector).
#' @param log_target A function taking numeric vector that evaluates the natural logarithm of the
#' unnormalized target density. Returns a numeric scalar.
#' @param mu A numeric vector with the location parameter of the Student t pseudo-target.
#' @param Sig Positive definite scale matrix of the Student t pseudo-target. Alternatively, a
#' lower-triangular matrix with the Cholesky factor of the scale matrix
#' (for faster computation).
#' @param df Degrees of freedom of Student t pseudo-target.
#'
#' @return A list contains two elements: \code{x} is the new state and \code{nEvaluations}
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @references
#' Nishihara, R., Murray, I., and Adams, R. P. (2014), "Parallel MCMC with Generalized Elliptical Slice Sampling," *Journal of Machine Learning Research*, 15, 2087-2112. \url{https://jmlr.org/papers/v15/nishihara14a.html}
#'
#' @importFrom stats dt rgamma
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' n_iter <- 10 # set to 1e4 for more complete illustration
#' draws <- matrix(0.3, nrow = n_iter, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, n_iter)) {
#'   out <- slice_genelliptical_mv(draws[i - 1,], log_target = lf,
#'               mu = c(0.5, 0.5), Sig = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2),
#'               df = 5)
#'   draws[i,] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (n_iter - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#'
slice_genelliptical_mv <- function(
  x,
  log_target,
  mu,
  Sig,
  df,
  is_chol = FALSE
) {
  K <- length(x)
  stopifnot(
    length(mu) == K,
    dim(Sig) == c(K, K),
    length(df) == 1L,
    all(is.finite(x)),
    is.function(log_target),
    all(is.finite(mu)),
    all(is.finite(Sig)),
    is.finite(df),
    df > 0.0,
    is.logical(is_chol)
  )

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
    is_chol <- TRUE
  }

  a <- 0.5 * (df + K)
  b <- 0.5 * (df + drop(crossprod(forwardsolve(SigL, (x - mu)))))
  s <- 1.0 / rgamma(1, shape = a, rate = b) # rate of gamma <=> shape of inv-gamma

  lff <- function(xx) {
    log_target(xx) +
      a * log1p(drop(crossprod(forwardsolve(SigL, (xx - mu)))) / df)
  }

  .slice_elliptical_mv_impl(
    x = x,
    log_lik = lff,
    mu = mu,
    Sig = sqrt(s) * SigL,
    is_chol = is_chol,
    sampler = "slice_genelliptical_mv"
  )
}
