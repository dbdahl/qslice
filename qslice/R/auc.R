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
#' u_samples <- rbeta(10e3, 4, 2)
#' auc(u = u_samples)
#' auc(u = u_samples, nbins = 50)
#' auc(y = function(x) {dbeta(x, 4, 2)}, nbins = 30)
#' auc(y = function(x) {dbeta(x, 4, 2)}, nbins = 300)
#' xx <- seq(0.001, 0.999, length = 1000)
#' auc(x = xx, y = function(x) {dbeta(x, 4, 2)})
#' auc(x = xx, y = dbeta(xx, 4, 2))
auc <- function(u = NULL, x = NULL, y = NULL, nbins = 30) {
  xchecks_auc <- function(x) {
    stopifnot(
      is.numeric(x),
      all(is.finite(x)),
      all(diff(x) > 0),
      x[1] > 0.0,
      x[length(x)] < 1.0
    )
  }

  if (is.null(u)) {
    if (is.function(y)) {
      if (is.null(x)) {
        x <- seq(1e-6, 1.0 - 1e-6, length = nbins) # reach into tails of function
      }

      y <- sapply(x, FUN = y) # convert y to density evaluations
      yn <- y / max(y)
      out <- mean(yn)
    } else {
      xchecks_auc(x)
      nx <- length(x)
      stopifnot(nx == length(y))

      breaks <- if (nx == 1L) {
        c(0, 1)
      } else {
        c(0, (x[-nx] + x[-1L]) / 2, 1)
      }

      widths <- diff(breaks)
      yn <- y / max(y)
      out <- sum(yn * widths)
    }
  } else {
    ## u supplied
    stopifnot(
      is.numeric(u),
      length(u) > 0L,
      all(is.finite(u)),
      all(u >= 0.0),
      all(u <= 1.0)
    )
    if (is.null(x)) {
      bins <- seq(0.0, 1.0, len = nbins + 1)
      x <- (bins[-(nbins + 1)] + bins[-1]) / 2.0
      y <- tabulate(as.numeric(cut(u, breaks = bins)), nbins = nbins)
      stopifnot(sum(y) > 0)
      yn <- y / max(y)
      out <- mean(yn)
    } else {
      xchecks_auc(x)
      nx <- length(x)
      breaks <- if (nx == 1L) {
        c(0, 1)
      } else {
        c(0, (x[-nx] + x[-1L]) / 2, 1)
      }
      binid <- findInterval(
        u,
        breaks,
        rightmost.closed = TRUE,
        all.inside = TRUE
      )
      y <- tabulate(binid, nbins = nx)
      y <- y / max(y)
      stopifnot(sum(y) > 0)
      out <- auc(x = x, y = y, nbins = nx) # recursively calculate
    }
  }
  out
}
