#' Utility for a given target and pseudo-target
#'
#' Takes a pseudo-target and target (or samples from the target) and
#' evaluates the utility function for the transformed target, which can be one of
#' Area Under the Curve (AUC) and Mean Slice Width (MSW). See Heiner et al. (2024+).
#'
#' Optionally plot the target and pseudo-target densities as well as the
#' transformed target.
#'
#' @param pseudo List containing the following functions relating to the pseudo-target, each with scalar input:
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
#' Use of "function" requires specification of \code{log_target}. Utility is calculated using numerical integration (\link[stats]{integrate}) and, additionally for AUC, optimization (\link[stats]{optim}).
#'
#' Use of "samples" requires specification of \code{samples}.
#'
#' Use of "grid" requires specification of \code{log_target} and optionally \code{x}. Grid method is useful if numerical optimization or integration fail employed with \code{type = "function"}.
#'
#' @param x Optional numeric vector specifying grid (on (0,1)) over which to evaluate
#' the transformed target. Defaults to \code{NULL}.
#' @param nbins Number of histogram bins to use (defaults to 30). Overridden by the length
#' of \code{x} if \code{x} is supplied.
#' @param plot Logical for whether to generate two plots:
#'
#' 1) direct comparison of the target and pseudo-target densities, and
#'
#' 2) transformed target density.
#'
#' Defaults to \code{TRUE}.
#' @param utility_type String identifying utility type, either AUC (default) or MSW.
#' @param options List containing the following. \code{tol_int}: Positive numeric scalar that passes to \code{abs.tol} in the call to \link[stats]{integrate} used for type = "function".
#' \code{n_opt}: Positive integer giving the number of equally spaced initialization points for numerical optimization (\link[stats]{optim}) in AUC with type = "function".
#' \code{interval_opt}: Numerical vector of length two giving the interval for numerical optimization (\link[stats]{optim}) in AUC with type = "function".
#'
#' @returns Scalar value of the utility function evaluation.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*.
#'
#' @importFrom graphics legend lines hist
#' @importFrom graphics points
#' @importFrom stats integrate optim
#' @export
#' @examples
#' pseu <- pseudo_list(family = "t", params = list(loc = 0.25, sc = 1.0, degf = 1))
#' ltarg <- list(ld = function(x) dnorm(x, log = TRUE))
#' oldpar <- par(mfrow = c(1,2))
#' utility_pseudo(pseudo = pseu, log_target = ltarg$ld, type = "function", utility_type = "MSW")
#' utility_pseudo(pseudo = pseu, log_target = ltarg$ld, type = "function", utility_type = "AUC")
#' samp <- rnorm(10e3)
#' utility_pseudo(pseudo = pseu, samples = samp, nbins = 30, type = "samples", utility_type = "MSW")
#' utility_pseudo(pseudo = pseu, samples = samp, nbins = 30, type = "samples", utility_type = "AUC")
#' par(oldpar)
#'
utility_pseudo <- function(
  pseudo,
  log_target = NULL,
  samples = NULL,
  type = "samples",
  x = NULL,
  nbins = 30,
  plot = TRUE,
  utility_type = "AUC",
  options = list(
    tol_int = 1.0e-3, # integration accuracy
    n_opt = 10, # number of inits for numerical optimization in AUC with type = "integration"
    interval_opt = c(0.000001, 0.999999) # interval for numerical optimization in AUC with type = "integration"
  )
) {
  if (type == "samples") {
    u <- vapply(
      seq_along(samples),
      function(i) {
        .eval_pseudo_cdf(
          p = pseudo$p,
          x = samples[i],
          sampler = "utility_pseudo",
          stage = "transforming target samples",
          attempt = i,
          expected_length = 1L,
          require_interior = TRUE
        )
      },
      numeric(1)
    )
  } else {
    h <- function(psi) exp(log_target(pseudo$q(psi)) - pseudo$ld(pseudo$q(psi)))
  }

  if (isTRUE(plot)) {
    x_plot <- seq(
      pseudo$params$loc - 4.0 * pseudo$params$sc,
      pseudo$params$loc + 4.0 * pseudo$params$sc,
      length = 200
    )
    y_pseu_plot <- sapply(x_plot, function(z) exp(pseudo$ld(z)))

    if (type %in% c("function", "grid")) {
      y_targ_plot <- sapply(x_plot, function(z) exp(log_target(z)))
      y_targ_plot <- y_targ_plot / max(y_targ_plot)
      y_pseu_plot <- y_pseu_plot / max(y_pseu_plot)

      plot(
        x_plot,
        y_targ_plot,
        type = "l",
        lwd = 2,
        xlab = "x",
        ylab = "density (unnormalized)",
        main = paste0("Pseudo-target:\n", pseudo$txt),
        cex.main = 0.97
      )
      lines(x_plot, y_pseu_plot, lwd = 2, col = "red", lty = 2)
      legend(
        "topleft",
        lwd = 2,
        lty = c(1, 2),
        col = c("black", "red"),
        bty = "n",
        legend = c("target", "pseudo-target")
      )
    } else if (type == "samples") {
      y_pseu_plot <- y_pseu_plot
      hist(
        samples,
        freq = FALSE,
        xlab = "Samples",
        main = paste0("Pseudo-target:\n", pseudo$txt),
        cex.main = 0.97
      )
      lines(x_plot, y_pseu_plot, lwd = 2, col = "red", lty = 2)
      legend(
        "topleft",
        lwd = 2,
        lty = c(1, 2),
        col = c("black", "red"),
        bty = "n",
        legend = c("target", "pseudo-target")
      )
    }
  }

  utility_shrinkslice(
    h = h,
    x = x,
    u = u,
    type = ifelse(type == "function", "integration", type),
    nbins = nbins,
    plot = plot,
    utility_type = utility_type,
    options = options
  )
}


#' Utility for a Transformed Target
#'
#' Evaluates the utility function for a transformed target, which can be one of
#' Area Under the Curve (AUC) and Mean Slice Width (MSW).
#'
#' @param h Function to evaluate the unnormalized transformed target
#' \eqn{h(\psi) = g(\hat{\Pi}^{-1}(\psi))/\hat{\pi}(\hat{\Pi}^{-1}(\psi))}
#' with argument \eqn{\psi \in (0,1)}. Use \code(h = NULL) if type = "samples".
#' @param x Numeric vector of histogram locations. Not used if \code{u} is supplied or type = "integration". If used but not supplied, a default grid is chosen.
#' @param u Numeric vector of samples supported on unit interval (\eqn{\psi}) with which to
#' create histogram (use \code{u = NULL} if type = "integration" or "grid").
#' @param type String specifying the utility calculation method. One of "samples", "integration", or "grid".
#' Use of "integration" requires \code{h}. Use of "samples" requires \code{u}.
#' Use of "grid" requires \code{h} and optionally \code{x}. Default is to use "samples".
#' @param nbins Number of histogram bins to use (defaults to 30). Overridden by the length of x.
#' @param plot Logical for whether to plot a visualization of the transformed target.
#' Defaults to \code{FALSE}.
#' @param utility_type String identifying utility type, either "AUC" (default) or "MSW"
#' @param options List containing the following. \code{tol_int}: Positive numeric scalar that passes to \code{abs.tol} in the call to \code{integrate()} used for type = "integration".
#' \code{n_opt}: Positive integer giving the number of initialization points for numerical optimization in AUC with type = "integration".
#' \code{interval_opt}: Numerical vector of length two giving the interval for numerical optimization in AUC with type = "integration".
#' @returns Scalar value of the utility function evaluation.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*.
#'
#' @importFrom graphics legend lines hist
#' @importFrom graphics points
#' @importFrom stats integrate optim
#' @export
#' @examples
#' a <- 4
#' b <- 2
#' u <- rbeta(10e3, a, b)
#' h <- function(x) dbeta(x, a, b)
#' x <- sort(runif(20))
#' type <- "AUC"
#' type <- "MSW"
#'
#' utility_shrinkslice(
#'   h = h,
#'   x = NULL,
#'   u = NULL,
#'   type = "integration",
#'   plot = TRUE,
#'   utility_type = type
#' )
#'
#' utility_shrinkslice(
#'   h = h,
#'   x = NULL,
#'   u = NULL,
#'   nbins = 50,
#'   type = "grid",
#'   plot = TRUE,
#'   utility_type = type
#' )
#'
#' utility_shrinkslice(
#'   h = h,
#'   x = x,
#'   u = NULL,
#'   type = "grid",
#'   plot = TRUE,
#'   utility_type = type
#' )
#'
#' utility_shrinkslice(
#'   h = NULL,
#'   x = NULL,
#'   u = u,
#'   nbins = 30,
#'   type = "samples",
#'   plot = TRUE,
#'   utility_type = type
#' )
#'
#' utility_shrinkslice(
#'   h = NULL,
#'   x = x,
#'   u = u,
#'   type = "samples",
#'   plot = TRUE,
#'   utility_type = type
#' )
#'
utility_shrinkslice <- function(
  h = NULL,
  x = NULL,
  u = NULL,
  type = "samples", # supplied with u. Alternatively, type = "samples", type = "integration" (supplied with function h) or "grid" (supplied with x, y)
  nbins = 30,
  plot = FALSE,
  utility_type = "AUC", # one of "AUC", "MSW"
  options = list(
    tol_int = 1.0e-3, # integration accuracy
    n_opt = 10, # number of inits for numerical optimization in AUC with type = "integration"
    interval_opt = c(0.000001, 0.999999) # interval for numerical optimization in AUC with type = "integration"
  )
) {
  if (type == "integration") {
    ## the supplied function here is the transformed target with support on (0, 1)

    stopifnot(is.function(h))
    if (isFALSE(is.null(x))) {
      warning(
        "Used type = 'integration' with utiliy_shrinkslice(). Performing numerical integration and discarding supplied x values."
      )
    }

    if (utility_type == "AUC") {
      util_out <- auc_int(
        h = h,
        tol_int = options$tol_int,
        n_opt = options$n_opt,
        interval_opt = options$interval_opt
      )
    } else if (utility_type == "MSW") {
      util_out <- meanSliceWidth_int(h = h, tol = options$tol_int)
    }

    if (isTRUE(plot)) {
      nbins <- 500
      x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins)
      y <- sapply(x, FUN = h)
    }
  } else {
    if (is.null(x)) {
      x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins) # important to reach into tails
    } else {
      nbins <- length(x)

      stopifnot(
        is.numeric(x),
        all(is.finite(x)),
        all(diff(x) > 0),
        x[1] > 0.0,
        x[length(x)] < 1.0
      )
    }

    breaks <- if (nbins == 1L) {
      c(0.0, 1.0)
    } else {
      c(0.0, (x[-nbins] + x[-1L]) / 2.0, 1.0)
    }
    widths <- diff(breaks)

    if (type == "samples") {
      stopifnot(
        is.numeric(u),
        length(u) > 0L,
        all(is.finite(u)),
        all(u >= 0.0),
        all(u <= 1.0)
      )
      binid <- findInterval(
        u,
        breaks,
        rightmost.closed = TRUE,
        all.inside = TRUE
      )

      counts <- tabulate(binid, nbins = nbins)
      y <- counts / widths # histogram height
      stopifnot(sum(y) > 0)
    } else if (type == "grid") {
      stopifnot(is.function(h))
      y <- sapply(x, FUN = h)
      stopifnot(all(y >= 0.0), any(y > 0.0))
    }

    if (utility_type == "AUC") {
      util_out <- sum(y * widths) / max(y) # auc
    } else if (utility_type == "MSW") {
      normc <- sum(y * widths)
      if (normc <= 0) {
        stop("MSW calculation: histogram has zero total mass.")
      }

      alpha <- vapply(
        y,
        function(hi) {
          sum(pmin(y, hi) * widths) / normc
        },
        numeric(1)
      )

      util_out <- sum(alpha * widths) # msw
    }
  }

  if (isTRUE(plot)) {
    y_plot <- y / max(y)

    plot(
      x,
      y_plot,
      xlim = c(-0.1, 1.1),
      ylim = c(0, max(y_plot)),
      type = ifelse(type == "integration", "l", "p"),
      lwd = 2,
      xlab = expression(psi),
      ylab = expression(h(psi)),
      main = paste0(
        "Transformed target\n",
        utility_type,
        ": ",
        round(util_out, 3)
      )
    )

    if (type != "integration") {
      lines(
        rep(breaks, each = 2),
        c(0.0, rep(y_plot, each = 2), 0.0),
        type = "s"
      )
    }
  }

  util_out
}
