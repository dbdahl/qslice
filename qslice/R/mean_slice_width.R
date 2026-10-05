# #' Expected Slice Width (integration)
# #'
# #' Calculate the expected slice width of a slice sampler for an unnormalized
# #' density supported on the unit interval.
# #'
# #' Uses numeric integration to evaluate a double integral.
# #'
# #'
# #' @param h Function evaluating an unnormalized density supported on the unit interval.
# #' @param tol Numerical scalar passed to \code{abs.tol} in the call to \code{integrate()}.
# #'
# #' @importFrom stats integrate
# #' @keywords internal
# #'
.meanSliceWidth_int <- function(h, tol = 0.005) {
  interval <- c(0.0, 1.0)

  h_vec <- function(x) sapply(x, FUN = h) # vectorized version of h

  nc <- integrate(h_vec, lower = interval[1], upper = interval[2])$value
  h_norm <- function(x) h_vec(x) / nc

  h_fill <- function(x, hnorm_at_outer) {
    # x must be able to be a vector
    pmin(h_norm(x), hnorm_at_outer)
  }

  inner_int <- function(x_outer) {
    y <- h_norm(x_outer)
    sapply(y, function(z) {
      # inner_int must be vectorized
      integrate(
        h_fill,
        lower = interval[1],
        upper = interval[2],
        hnorm_at_outer = z,
        abs.tol = tol
      )$value
    })
  }

  outer_int <- tryCatch(
    {
      integrate(
        inner_int,
        lower = interval[1],
        upper = interval[2],
        abs.tol = tol
      )
    },
    error = function(e) {
      message(
        "\nUnable to calculate MSW with numerical integration.\n
        If optimizing the pseudo-target, try different initial parameter values.\n
        Consider calculating MSW on a grid (use type = 'grid') instead.\n
        Error message from integrate():\n"
      )

      stop(e)
    }
  )

  outer_int$value
}
