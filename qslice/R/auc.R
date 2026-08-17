# #' Area Under the Curve (integration)
# #'
# #' Calculate the area under the curve after restricting
# #' the curve to fit within the unit square. Specifically, the highest unnormalized density value reaches 1 and
# #' the support is the unit interval. See Heiner et al. (2024+).
# #'
# #' Uses optimization (\code{optim()}) to find the highest point and numeric integration (\code{integrate()}) to find the area.
# #'
# #'
# #' @param h Function evaluating an unnormalized density supported on the unit interval.
# #' @param tol_int Numerical scalar passed to \code{abs.tol} in the call to \code{integrate()}.
# #' @param n_opt Positive integer specifying the number of initializtion points (equally spaced over (0,1)) for the optimization routine.
# #' @param interval_opt Numeric vector of length two (in (0,1)) giving the interval over which to optimize. Defaults to \code{c(0.000001, 0.999999)}.
# #'
# #' @importFrom stats integrate
# #' @keywords internal
# #'
auc_int <- function(
  h,
  tol_int = 0.005,
  n_opt = 10,
  interval_opt = c(0.000001, 0.999999)
) {
  stopifnot(is.function(h))
  stopifnot(n_opt > 0)
  stopifnot(tol_int >= 0.0)

  inits <- seq(0.01, 0.99, length = n_opt)
  stopifnot(
    interval_opt[1] < min(inits),
    interval_opt[2] > max(inits),
    interval_opt[1] >= 0.0,
    interval_opt[2] <= 1.0,
    interval_opt[1] < interval_opt[2]
  )

  opt <- lapply(inits, function(x) {
    result <- tryCatch(
      {
        optim(
          x,
          h,
          control = list(fnscale = -1),
          lower = interval_opt[1],
          upper = interval_opt[2],
          method = "L-BFGS-B"
        )
      },
      error = function(e) {
        out <- list()
        out$convergence <- 100 # anything but 0
        out$message <- e$message
        out
      }
    )
  })

  convergence <- sapply(opt, function(x) x$convergence)

  if (all(convergence != 0)) {
    stop(
      "AUC optimiation for unnormalized target density failed to converge.\n 
      Consider calculating AUC on a grid (use type = 'grid') instead of via numerical optimization/integration.\n
      Error messages from optimization step:\n",
      paste(sapply(opt, function(x) x$message), sep = "\n")
    )
  } else {
    indx_converge <- which(convergence == 0)
    opt_vals <- sapply(indx_converge, function(j) opt[[j]]$value)
    max_val <- max(opt_vals)
  }

  h_vec <- function(x) sapply(x, FUN = h) # vectorized version of h

  int <- integrate(
    h_vec,
    lower = 0.0,
    upper = 1.0,
    abs.tol = tol_int
  )

  auc_val <- int$value / max_val

  if (auc_val > (1.0 + int$abs.error)) {
    stop(
      "Calculated AUC is greater than 1.\nAUC = ",
      auc_val,
      ", integral = ",
      int$value,
      ", max(h(x)) = ",
      max_val,
      "\nConsider calculating AUC on a grid (use type = 'grid') instead of via numerical optimization/integration."
    )
  } else if (auc_val < (0.0 - int$abs.error)) {
    stop(
      "Calculated AUC is less than 0.  AUC = ",
      auc_val,
      ", integral = ",
      int$value,
      ", max(h(x)) = ",
      max_val,
      "\nConsider calculating AUC on a grid (use type = 'grid') instead of via numerical optimization/integration."
    )
  }

  min(auc_val, 1.0)
}
