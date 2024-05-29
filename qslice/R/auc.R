#' Area Under the Curve (histogram)
#'
#' Calculate the histogram approximation to the area under the curve after restricting
#' the curve to fit within the unit square. Specifically, the highest histogram bar reaches 1 and
#' the support is the unit interval. See Heiner et al. (2024+).
#'
#' Accepts either samples \code{u} or a function \code{y} representing a (possibly
#' unnormalized) probability density supported on the unit interval.
#'
#' @param u Numeric vector of samples supported on unit interval with which to
#' create a histogram (use \code{u = NULL} if \code{x} and \code{y} are supplied).
#' @param x Numeric vector of histogram locations. (Not used if \code{u} is supplied).
#' @param y Numeric vector of histogram heights OR function evaluating the curve
#' for a given value of \code{u} supported on (0,1). (Not used if \code{u} is supplied).
#' @param nbins Number of histogram bins to use (defaults to 30).
#'
#' @returns The (approximate) area under the curve as a numeric value of length one.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*.
#'
#' @export
#' @examples
#' u_samples <- rbeta(10e3, 2, 2)
#' auc(u = u_samples)
#' auc(u = u_samples, nbins = 50)
#' auc(y = function(x) {dbeta(x, 2, 2)}, nbins = 30)
#' auc(y = function(x) {dbeta(x, 2, 2)}, nbins = 300)
#' xx <- seq(0.001, 0.999, length = 1000)
#' auc(x = xx, y = function(x) {dbeta(x, 2, 2)})
#' auc(x = xx, y = dbeta(xx, 2, 2))
auc <- function(u = NULL, x = NULL, y = NULL, nbins = 30) {

  if (is.null(u)) {

    if (is.function(y)) {

      if (is.null(x)) {
        x <- seq(1e-6, 1.0 - 1e-6, length = nbins)
      }

      y <- sapply(x, FUN = y)

    }

    stopifnot(length(x) == length(y))
    stopifnot(all(x > 0.0) && all(y >= 0.0) && all(x < 1.0))

  } else {

    bins <- seq(0.0, 1.0, len = nbins + 1)
    x <- (bins[-(nbins+1)] + bins[-1]) / 2.0
    y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins = nbins)

  }

  yn <- y / max(y)

  mean(yn)
}
