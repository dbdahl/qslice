#' Calculate the Utility Function for a Transformed Target
#'
#' Evaluates the utility function for a transformed target. The base utility can be one of
#' Area Under the Curve (AUC) and Expected Slice Width. Optionally adds a penalty for multimodality,
#' measured as the area of water the density could hold.
#'
#' @param h Function to evaluate the unnormalized transformed target \eqn{h(\psi) = g(\hat{\Pi}^{-1}(\psi))/\hat{\pi}(\hat{\Pi}^{-1}(\psi))} with argument \eqn{\psi \in (0,1)}.
#' @param x Numeric vector of histogram locations. (Not used if \code{u} is supplied).
#' @param y Numeric vector of histogram heights OR function evaluating the curve
#' for a given value of \code{u}.(Not used if \code{u} is supplied).
#' @param u Numeric vector of samples supported on unit interval (\eqn{\psi}) with which to
#' create histogram (use \code{u = NULL} if \code{x} and \code{y} are supplied).
#' @param type String specifying the input type. One of "function", "samples", or "grid".
#' Use of "function" requires \code{h}. Use of "samples" requires \code{samples}.
#' Use of "grid" requires \code{h}. Defaults to "samples".
#' @param coeffs Positive numeric vector of length two giving relative weights of
#' 1-the base utility (one of AUC or mean slice width) and 2-the penalty for multimodality,
#' measured as the area of water the density could hold. Defaults to \code{c(1.0, 0.0)}.
#' @param nbins Number of histogram bins to use (defaults to 30).
#' @param plot Logical for whether to plot a visualization of the transformed target.
#' Defaults to \code{TRUE}.
#' @param use_meanSliceWidth Logical for whether the base utility should use
#' Expected Slice Width (if \code{TRUE}; slower) or Area Under the Curve (AUC) if
#' \code{FALSE}. Defaults to \code{FALSE}.
#' @param tol_int Positive numeric scalar that passes to \code{abs.tol} in the call to \code{integrate()}.
#' Defaults to \code{1.0e-3}.
#' @returns Named vector with:
#'
#'  \code{util}: final utility calculated as \code{coeffs[1]*base_util - coeffs[2]*water_area}
#'
#'  \code{auc}: area of the unit square that is under the curve \eqn{h}
#'
#'  \code{water_area}: area of the unit square that could hold water (convex portions of \eqn{h})
#'
#'  \code{msw}: Mean slice width (if used as the base utility)
#'
#' @export
#' @examples
#' x <- seq(0.001, 0.999, length = 1000)
#' y <- dbeta(x, 2.0, 2.0)
#' utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE, use_meanSliceWidth = TRUE)
#' utility_shrinkslice(h = \(x) dbeta(x, 2.0, 2.0), type = "function", plot = TRUE, use_meanSliceWidth = TRUE)
#' utility_shrinkslice(u = rbeta(1e3, 2.0, 2.0), type = "samples")
utility_shrinkslice <- function(h = NULL, x = NULL, y = NULL, u = NULL,
                                type = "samples", # supplied with u. Alternatively, type = "samples_kde", type = "function" (supplied with function h) or "grid" (supplied with x, y)
                                coeffs = c(1.0, 0.0), # default is to have no water penalty
                                nbins = 30,
                                plot = FALSE,
                                use_meanSliceWidth = FALSE,
                                tol_int = 1.0e-3) {

  if (type %in% c("function", "samples_kde")) {

    ## the supplied function here is the transformed target with support on (0, 1)

    if (type == "samples_kde") {
      h = beta_kde(u)
    }

    tmp <- water_area_int(h, interval = c(0.0 + 1.0e-9, 1.0 - 1.0e-9),
                          plot = plot, eps = 1.0e-3, tol_int = tol_int,
                          title = TRUE)
    auc <- tmp$AUC / tmp$totalArea
    wtr <- tmp$totalWaterArea / tmp$totalArea

    if (use_meanSliceWidth) {
      msw <- meanSliceWidth_int(h, tol = tol_int)
    }

  } else if (type %in% c("grid", "samples")) {

    if (type == "samples") {

      if (is.null(x)) {
        x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins) # trouble if it doesn't reach far enough into tails
      } else {
        stopifnot(length(x) == nbins)
      }

      bins <- c(0.0, x[-c(1, nbins)], 1.0)
      y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins = nbins )

      stopifnot(sum(y) > 0)

    }

    auc <- auc(x = x, y = y)
    wtr <- water_area(x = x, y = y, plot = plot)

    if (use_meanSliceWidth) {
      msw <- meanSliceWidth_grid(x = x, y = y)
    }

  } # currently no support for histogram method (unless supplied as x and y)

  if (use_meanSliceWidth) {
    out <- c( util = coeffs %*% c(msw, -wtr), auc = auc, water_area = wtr, msw = msw )
  } else {
    out <- c( util = coeffs %*% c(auc, -wtr), auc = auc, water_area = wtr )
  }

  out
}

util_pseu <- function(pseu, target = NULL, samples = NULL,
                      type = "samples",
                      x = NULL,
                      bins = NULL, nbins = NULL,
                      coeffs = c(1.0, 0.0),
                      plot = FALSE,
                      use_meanSliceWidth = TRUE,
                      tol_int = 1.0e-3) {

  if (type == "function") {

    h <- function(psi) exp( target$ld( pseu$q(psi) ) - pseu$ld( pseu$q(psi) ) )
    x <- NULL
    y <- NULL
    u <- NULL

  } else {

    h <- NULL

    if (type == "grid") {

      xx <- pseu$q(x)
      ly <- target$ld( xx ) - pseu$ld( xx )
      y <- exp(ly - max(ly))
      u <- NULL

    } else if (type == "samples") {

      u <- pseu$p(samples)
      y <- NULL

    } else if (type == "samples_kde") {

      u <- pseu$p(samples)
      x <- NULL
      y <- NULL

    }

  }

  utility_shrinkslice(h = h, x = x, y = y, u = u,
                      type = type,
                      coeffs = coeffs,
                      nbins = nbins,
                      plot = plot,
                      use_meanSliceWidth = use_meanSliceWidth,
                      tol_int = tol_int)
}
