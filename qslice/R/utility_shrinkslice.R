#' Utility for a given target and pseudo-target
#'
#' Takes a pseudo-target and target (or samples from the target) and
#' evaluates the utility function for the transformed target, which can be one of
#' Area Under the Curve (AUC) and Mean Slice Width (MSW). See Heiner et al. (2024+).
#'
#' Optionally plot the target and pseudo-target densities as well as the
#' transformed tartet.
#'
#' @param pseudo List containing the following functions with scalar input:
#'
#'  \code{ld}: function to evaluate the log density
#'
#'  \code{q}: function to evaluate the quantile function
#'
#'  \code{p}: function to evaluate the distribution function
#'
#' @param log_target Function to evaluate the log density of the unnormalized target.
#'
#' @param samples Numeric vector of samples from the target distribution.
#' @param type String specifying the input type. One of "function", "samples", or "grid".
#' Default is to use "samples".
#'
#' Use of "function" requires specification of \code{log_target}.
#'
#' Use of "samples" requires specification of \code{samples}.
#'
#' Use of "grid" requires specification of \code{x}.
#'
#' @param x Numeric vector specifying grid (on (0,1)) over which to evaluate
#' the transformed target. Defaults to \code{NULL}.
#' @param nbins Number of histogram bins to use (defaults to 30). Must match the length
#' of \code{x} if \code{x} is supplied.
#' @param plot Logical for whether to generate two plots:
#'
#' 1) direct comparison of the target and pseudo-target densities, and
#'
#' 2) transformed target density.
#'
#' Defaults to \code{TRUE}.
#' @param utility_type String identifying utility type, either AUC (default) or MSW.
#' @param tol_int Positive numeric scalar that passes to \code{abs.tol} in the call to \link[stats]{integrate}.
#' Defaults to \code{1.0e-3}.
#' @returns Scalar value of the utility function evaluation.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*.
#'
#' @importFrom graphics legend lines
#' @importFrom graphics curve points segments text
#' @importFrom stats integrate
#' @export
#' @examples
#' pseu <- pseudo_list(family = "logistic", params = list(loc = 0.0, sc = 0.66))
#' ltarg <- list(ld = function(x) dnorm(x, log = TRUE))
#' oldpar <- par(mfrow = c(1,2))
#' utility_pseudo(pseudo = pseu, log_target = ltarg$ld, type = "function",
#'                nbins = 100, utility_type = "MSW")
#' samp <- rnorm(10e3)
#' utility_pseudo(pseudo = pseu, samples = samp, type = "samples", utility_type = "AUC")
#' utility_pseudo(pseudo = pseu, samples = samp, type = "samples", utility_type = "MSW")
#' par(oldpar)
#'
utility_pseudo <- function(pseudo, log_target = NULL, samples = NULL,
                        type = "samples",
                        x = NULL,
                        nbins = 30,
                        plot = TRUE,
                        utility_type = "AUC",
                        tol_int = 1.0e-3) {

  if (type == "function") {

    h <- function(psi) exp( log_target( pseudo$q(psi) ) - pseudo$ld( pseudo$q(psi) ) )
    x <- NULL
    y <- NULL
    u <- NULL

  } else {

    h <- NULL

    if (type == "grid") {

      xx <- pseudo$q(x)
      ly <- sapply(xx, FUN = function(z) log_target( z ) - pseudo$ld( z ))
      y <- exp(ly - max(ly))
      u <- NULL

    } else if (type == "samples") {

      u <- sapply(samples, FUN = pseudo$p)
      y <- NULL

    }

  }

  if (isTRUE(plot)) {

    if (type == "function") {
      x_plot <- seq(pseudo$params$loc - 3.5*pseudo$params$sc,
                    pseudo$params$loc + 3.5*pseudo$params$sc,
                    length = 100)

      y_targ_plot <- sapply(x_plot, function(z) exp(log_target(z)))
      y_targ_plot <- y_targ_plot / max(y_targ_plot)
      y_pseu_plot <- sapply(x_plot, function(z) exp(pseudo$ld(z)))
      y_pseu_plot <- y_pseu_plot / max(y_pseu_plot)

      plot(x_plot, y_targ_plot, type = "l", lwd = 2,
           xlab = "x", ylab = "density (unnormalized)",
           main = paste0("Pseudo-target:\n", pseudo$txt))
      lines(x_plot, y_pseu_plot, lwd = 2, col = "red")
      legend("topleft", lwd = 2, col = c("black", "red"), bty = "n",
             legend = c("target", "pseudo-target"))
    }
  }

  utility_shrinkslice(h = h, x = x, y = y, u = u,
                      type = type,
                      nbins = nbins,
                      plot = plot,
                      utility_type = utility_type,
                      tol_int = tol_int)
}



# #' Utility for a Transformed Target
# #'
# #' Evaluates the utility function for a transformed target, which can be one of
# #' Area Under the Curve (AUC) and Mean Slice Width (MSW).
# #'
# #' @param h Function to evaluate the unnormalized transformed target
# #' \eqn{h(\psi) = g(\hat{\Pi}^{-1}(\psi))/\hat{\pi}(\hat{\Pi}^{-1}(\psi))}
# #' with argument \eqn{\psi \in (0,1)}.
# #' @param x Numeric vector of histogram locations. (Not used if \code{u} is supplied).
# #' @param y Numeric vector of histogram heights OR function evaluating the curve
# #' for a given value of \code{u}.(Not used if \code{u} is supplied).
# #' @param u Numeric vector of samples supported on unit interval (\eqn{\psi}) with which to
# #' create histogram (use \code{u = NULL} if \code{x} and \code{y} are supplied).
# #' @param type String specifying the input type. One of "function", "samples", or "grid".
# #' Use of "function" requires \code{h}. Use of "samples" requires \code{u}.
# #' Use of "grid" requires \code{h} and optionally \code{x}. Default is to use "samples".
# #' @param nbins Number of histogram bins to use (defaults to 30).
# #' @param plot Logical for whether to plot a visualization of the transformed target.
# #' Defaults to \code{FALSE}.
# #' @param utility_type String identifying utility type, either AUC (default) or MSW
# #' @param tol_int Positive numeric scalar that passes to \code{abs.tol} in the call to \code{integrate()}.
# #' Defaults to \code{1.0e-3}.
# #' @returns Scalar value of the utility function evaluation.
# #' @keywords internal
# #'
utility_shrinkslice <- function(h = NULL, x = NULL, y = NULL, u = NULL,
                                type = "samples", # supplied with u. Alternatively, type = "samples", type = "function" (supplied with function h) or "grid" (supplied with x, y)
                                nbins = 30,
                                plot = FALSE,
                                utility_type = "AUC",
                                tol_int = 1.0e-3) {

  if (type == "function") {

    ## the supplied function here is the transformed target with support on (0, 1)

    if (utility_type == "AUC") {

      util_out <- auc(y = h, nbins = nbins)

    } else if (utility_type == "MSW") {

      util_out <- meanSliceWidth_int(h = h, tol = tol_int)

    }

    if (isTRUE(plot)) {

      x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins) # trouble if it doesn't reach far enough into tails
      y <- sapply(x, FUN = h)
      y <- y / max(y)

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
      y <- y / max(y)

      stopifnot(sum(y) > 0)

    }

    if (utility_type == "AUC") {

      util_out <- auc(x = x, y = y)

    } else if (utility_type == "MSW") {

      util_out <- meanSliceWidth_grid(x = x, y = y)

    }

  } # currently no support for histogram method (unless supplied as x and y)

  if (isTRUE(plot)) {
    plot(x, y, type = "l", lwd = 2,
         xlab = expression(psi), ylab = expression(h(psi)),
         main = paste0("Transformed target\n", utility_type, ": ",
                       round(util_out,3))
         )
  }

  util_out
}
