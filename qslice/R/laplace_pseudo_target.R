#' Pseudo-target from Laplace Approximation
#'
#' Find the location and scale for an approximating pseudo-target via Laplace
#' approximation.
#'
#' @param log_target Univariate function evaluating the unnormalized log density to approximate.
#' @param init Numeric scalar for an initial value (used in optimization).
#' @param family String specifying the family of distributions for the pseudo-target.
#' Can be any of the families accepted by \link[qslice]{pseudo_list}.
#' @param params List specifying the parameters for the pseudo-target to be used.
#' The location and scale parameters will be replaced with the Laplace approximation and
#' others (e.g., degrees of freedom) will be retained.
#' @param sc_adj Positive numeric scalar; manual multiplicative adjustment to the
#' scale of the output pseudo-target.
#' @param lb Numeric scalar giving the value of left truncation of the resulting pseudo-target.
#' Defaults to \code{-Inf}.
#' @param ub Numeric scalar giving the value of right truncation of the resulting pseudo-target.
#' Defaults to \code{Inf}.
#' @param maxit See \link[stats]{optim}.
#' @param ... See \link[stats]{optim}.
#' @returns A list with the same outputs as \link[qslice]{pseudo_list}; also includes
#' \code{opt}, which gives output of \link[stats]{optim}.
#'
#' @importFrom stats optim
#'
#' @export
#' @examples
#' pseu <- lapprox(function(x) dnorm(x, log = TRUE),
#'   family = "t",
#'   params = list(loc = NA, sc = NA, degf = 5.0),
#'   init = 0.5, lb = -1.0)
#' curve(dnorm(x)/(1- pnorm(-1)), from = -1, to = 6, col = "blue")
#' xx <- seq(-1, 6, length = 500)
#' lines(xx, sapply(xx, FUN = pseu$d))
lapprox <- function(log_target, init, family = "t", params = NULL, sc_adj = 1.0,
                    lb = -Inf, ub = Inf, maxit = 100, ...) {

  fit <- optim(par = init, fn = log_target, control = list(fnscale = -1, maxit = maxit), method = 'BFGS', ...)
  loc <- fit$par
  hessian <- second_derivative( x = loc, h = 1e-5, f = log_target, ...)

  if (hessian > 0.0) warning(paste0("Hessian =", hessian, "; optim iters: ", fit$counts))

  sc <- sc_adj / sqrt(-hessian)

  params$loc <- loc
  params$sc <- sc

  out <- pseudo_list(family = family, params = params,
                     lb = lb, ub = ub, name = 'Laplace')

  out$opt <- fit

  out
}

second_derivative <- function(x, h = 1e-5, f, ... ) {

  num <- f(x + h, ...) - 2*f(x, ...) + f(x - h, ...)
  denom <- h^2

  num/denom
}
