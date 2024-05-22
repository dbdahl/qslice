
#' Multivariate Elliptical Slice Sampler
#'
#' Algorithm 1 of Nishihara et al. (2014) of the
#' elliptical slice sampler of Murray et al. (2010).
#'
#' @inherit slice_stepping_out
#' @param mu Numeric vector with the mean of the supporting normal distribution.
#' @param Sig Positive definite covariance matrix. Alternatively, a
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
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' n_iter <- 10 # set to 1e3 for more complete illustration
#' draws <- matrix(0.3, nrow = n_iter, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, n_iter)) {
#'   out <- slice_elliptical_mv(draws[i - 1,], log_target = lf,
#'               mu = c(0.5, 0.5), Sig = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2))
#'   draws[i,] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / (n_iter - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#'
slice_elliptical_mv <- function(x, log_target, mu, Sig, is_chol = FALSE) {
  nEvaluations <- 0

  k <- length(x)
  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x)
  }
  # Step 1
  fx <- f(x)
  stopifnot(fx > -Inf)
  y <- log(runif(1)) + fx

  nu <- SigL %*% rnorm(k, 0.0, 1.0) + mu

  twopi <- 2.0 * pi
  theta <- runif(1, 0, twopi)
  theta_min <- theta - twopi
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

#' Generalized Elliptical Slice Sampler (Multivariate)
#'
#' Generalized Elliptical Slice Sampler, Algorithm 2 of Nishihara et al. (2014)
#'
#' @inheritParams slice_elliptical_mv
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
slice_genelliptical_mv <- function(x, log_target, mu, Sig, df, is_chol = FALSE) {

  k <- length(x)
  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  a <- 0.5 * (df + k)
  b <- 0.5 * (df + drop(crossprod(forwardsolve(SigL, (x - mu)))))
  s <- 1.0 / rgamma(1, shape = a, rate = b) # rate of gamma <=> shape of inv-gamma

  lff <- function(xx) {
    log_target(xx) + a*log1p(drop(crossprod(forwardsolve(SigL, (xx - mu))))/df)
  }

  slice_elliptical_mv(x = x, log_target = lff, mu = mu, Sig = sqrt(s) * SigL,
                      is_chol = is_chol)
}
