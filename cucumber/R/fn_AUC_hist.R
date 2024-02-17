#' Area Under the Curve (histogram)
#'
#' Calculate the histogram approximation to the area under the curve after restricting
#' the curve to fit within the unit square. Specifically, the highest histogram bar reaches 1 and
#' the support is the unit interval.
#'
#'
#' @param x Numeric vector of histogram locations. (Not used if \code{u} is supplied).
#' @param y Numeric vector of histogram heights OR function evaluating the curve
#' for a given value of \code{u}.(Not used if \code{u} is supplied).
#' @param u Numeric vector of samples supported on unit interval with which to
#' create histogram (use \code{u = NULL} if \code{x} and \code{y} are supplied).
#' @param nbins Number of histogram bins to use (defaults to 30).
#'
#' @export
#' @examples
#' auc(u = rbeta(1000, 2, 2))
#' auc(x = runif(1000), y = function(x) {dbeta(x, 2, 2)})
auc <- function(x, y, u = NULL, nbins = 30) {

  if (is.null(u)) {

    if (is.function(y)) {
      y <- y(x)
    }

    stopifnot(length(x) == length(y))
    stopifnot(all(x > 0.0) && all(y >= 0.0) && all(x < 1.0))

  } else {

    bins = seq(0.0, 1.0, len = nbins + 1)
    x <- (bins[-(nbins+1)] + bins[-1]) / 2.0
    y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins = nbins)

  }

  yn <- y / max(y)

  mean(yn)
}

