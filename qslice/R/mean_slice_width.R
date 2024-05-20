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
meanSliceWidth_int <- function(h, tol = 0.005) {

  interval <- c(0.0, 1.0)

  h_vec <- function(x) sapply(x, FUN = h)

  nc <- integrate(h_vec, lower = interval[1], upper = interval[2])$value
  h_norm <- function(x) h_vec(x) / nc

  h_fill <- function(x, hnorm_at_outer) { # x must be able to be a vector
    pmin(h_norm(x), hnorm_at_outer)
  }

  inner_int <- function(x_outer) {
    y <- h_norm(x_outer)
    sapply(y, function(z) { # inner_int must be vectorized
      integrate(h_fill, lower = interval[1], upper = interval[2], hnorm_at_outer = z, abs.tol = tol)$value
    })
  }

  integrate(inner_int, lower = interval[1], upper = interval[2], abs.tol = tol)$value
}


# #' Expected Slice Width (over a grid)
# #'
# #' Calculate the expected slice width of a slice sampler for an unnormalized
# #' density supported on the unit interval.
# #'
# #' Uses an approximation to the density on a grid of \code{(x, y)} pairs.
# #'
# #'
# #' @param x Numeric vector of ordered values ranging from 0 to 1.
# #' @param y Numeric vector of evaluations of the unnormalized density at the
# #' supplied values of \code{x}.
# #'
# #' @keywords internal
# #'
meanSliceWidth_grid <- function(x, y) {
  n <- length(x)

  stopifnot(length(y) == n)

  sumy <- sum(y)
  yn <- y / sum(y)

  ord <- order(y, y + (1.0e-9)*runif(n, min = -1.0, max = 1.0)) # to break ties... not necessary?
  yn_ord <- yn[ord]
  A_inYorder <- (n:1)*yn_ord + c(0, cumsum(yn_ord[-n]))

  mean(A_inYorder)
}
