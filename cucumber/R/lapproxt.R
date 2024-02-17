
#' Helper Function for Laplace Approximation
#'
#'
second_derivative <- function( x, h = 1e-5, f ) {

  num <- f(x + h) - 2*f(x) + f(x - h)
  denom <- h^2

  num/denom
}

#' Cauchy Pseudo-Target from Laplace Approximation
#'
#' Author: Sam Johnson
#'
#' @param f Univariate function evaluating the unnormalized density to approximate.
#' Note that this is not on the log scale.
#' @param init Numeric scalar for an initial value (used in optimization).
#' @param sc_adj Positive numeric scalar; manual multiplicative adjustment to the
#' scale of the output pseudo-target.
#' @param lb Numeric scalar giving the value of left truncation of the resulting pseudo-target.
#' Defaults to \code{-Inf}.
#' @param ub Numeric scalar giving the value of right truncation of the resulting pseudo-target.
#' Defaults to \code{Inf}.
#'
#' @export
#' @examples
#' pseu <- lapproxt(function(x) exp(-x^2), init = 0.5, lb = 0.0)
#' curve(pseu$d(x), from = -2, to = 6)
lapproxt <- function(f, init, sc_adj = 1.0, lb = -Inf, ub = Inf, ...) {

  fit <- optim(par = init, fn = f, control = list(fnscale = -1), method = 'BFGS')
  loc <- fit$par
  hessian <- second_derivative( x = loc, h = 1e-5, f = f )
  sc <- sc_adj / sqrt(-hessian)
  out <- pseudo_t_list(loc = loc, sc = sc, degf = 1, lb = lb, ub = ub, name = 'Laplace')

  out
}
