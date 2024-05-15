
#' Multivariate Elliptical Slice Sampler
#'
#' Algorithm 1 of Nishihara et al (2014) of the
#' elliptical slice sampler of Murray, Adams, MacKay (2010).
#'
#' @inherit slice_stepping_out
#' @param mu Numeric vector with the mean of the normal to be sampled.
#' @param Sig Positive definite covariance matrix. Alternatively, a
#' lower-triangular matrix with the Cholesky factor of the covariance matrix
#' (for faster computation).
#' @param is_chol Logical, is the supplied \code{Sig} in Cholesky (lower triangular) format?
#'
#' @importFrom stats runif rnorm
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' draws <- matrix(0.3, nrow = 1000, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, nrow(draws))) {
#'   out <- slice_elliptical_mv(draws[i - 1,], target = lf,
#'               mu = c(0.5, 0.5), Sig = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2))
#'   draws[i,] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / nrow(draws)
#' nEvaluations / ess(draws[,1])
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#'
slice_elliptical_mv <- function(x, target, mu, Sig, is_chol = FALSE) {
  nEvaluations <- 0

  k <- length(x)
  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  # is_chol <- all(Sig[upper.tri(Sig, diag = FALSE)] == 0)
  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
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

#' General Elliptical Slice Sampler (Multivariate)
#'
#' General Elliptical Slice Sampler, Algorithm 2 of Nishihara (2014)
#'
#' @inheritParams slice_elliptical_mv
#' @param df Degrees of freedom of Student t pseudo-target.
#'
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @importFrom stats dt rgamma
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' draws <- matrix(0.3, nrow = 10e3, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, nrow(draws))) {
#'   out <- slice_genelliptical_mv(draws[i - 1,], target = lf,
#'               mu = c(0.5, 0.5), Sig = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2),
#'               df = 5)
#'   draws[i,] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#' }
#' nEvaluations / nrow(draws)
#' nEvaluations / ess(draws[,1])
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#'
slice_genelliptical_mv <- function(x, target, mu, Sig, df, is_chol = FALSE) {

  k <- length(x)
  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  # is_chol <- all(Sig[upper.tri(Sig, diag = FALSE)] == 0)
  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  a <- 0.5 * (df + k)
  b <- 0.5 * (df + drop(crossprod(forwardsolve(SigL, (x - mu)))))
  s <- 1.0 / rgamma(1, shape = a, rate = b) # rate of gamma <=> shape of inv-gamma

  lff <- function(xx) {
    target(xx) + a*log1p(drop(crossprod(forwardsolve(SigL, (xx - mu))))/df)
  }

  slice_elliptical_mv(x = x, target = lff, mu = mu, Sig = sqrt(s) * SigL,
                      is_chol = is_chol)
}