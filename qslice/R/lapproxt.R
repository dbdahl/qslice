#' Cauchy Pseudo-Target from Laplace Approximation
#'
#' Author: Sam Johnson
#'
#' @param lf Univariate function evaluating the unnormalized log density to approximate.
#' @param init Numeric scalar for an initial value (used in optimization).
#' @param sc_adj Positive numeric scalar; manual multiplicative adjustment to the
#' scale of the output pseudo-target.
#' @param lb Numeric scalar giving the value of left truncation of the resulting pseudo-target.
#' Defaults to \code{-Inf}.
#' @param ub Numeric scalar giving the value of right truncation of the resulting pseudo-target.
#' Defaults to \code{Inf}.
#' @param maxit See \link[coda]{optim}.
#' @param ... See \link[coda]{optim}.
#'
#' @export
#' @importFrom stats optim
#' @examples
#' pseu <- lapproxt(function(x) dnorm(x, log = TRUE), init = 0.5, lb = -1.0)
#' curve(dnorm(x)/(1- pnorm(-1)), from = -1, to = 6, col = "blue")
#' curve(pseu$d(x), from = -2, to = 6, add = TRUE)
lapproxt <- function(lf, init, sc_adj = 1.0, lb = -Inf, ub = Inf, maxit = 100, ...) {

  fit <- optim(par = init, fn = lf, control = list(fnscale = -1, maxit = maxit), method = 'BFGS', ...)
  loc <- fit$par
  hessian <- second_derivative( x = loc, h = 1e-5, f = lf, ...)

  if (hessian > 0.0) cat("Hessian =", hessian, "; optim iters: ", fit$counts, "\n")

  sc <- sc_adj / sqrt(-hessian)
  out <- pseudo_t_list(loc = loc, sc = sc, degf = 1, lb = lb, ub = ub, name = 'Laplace')

  out
}

second_derivative <- function( x, h = 1e-5, f, ... ) {

  num <- f(x + h, ...) - 2*f(x, ...) + f(x - h, ...)
  denom <- h^2

  num/denom
}
