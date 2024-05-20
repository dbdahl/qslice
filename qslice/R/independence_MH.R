#' Independence Metropolis-Hastings
#'
#' @param x The current state (scalar or numeric vector).
#' @param log_target A function taking a scalar or numeric vector that evaluates the log-target
#'   density, returning a numeric scalar.
#' @param pseudo List specifying the pseudo-target (proposal distribution). If the list length is
#' equal to the number of dimensions in \code{x}, each element is itself a list that specifies
#' the pseudo-target for the corresponding dimension with functions \code{ld}
#' that evaluates the log density for that dimension,
#' and \code{q} that evaluates the quantile (inverse-CDF) function for that dimension.
#' If the dimension of \code{x} is one, then supply only the inner list
#' specifying the single pseudo-target.
#'
#' If \code{x} is a vector but a single pseudo-target is supplied, the list must
#' contain a log-density function \code{ld} that accepts a vector, and a \code{r}
#' function that takes no arguments and generates a single multivariate draw from the
#' proposal distribution.
#' @return A list containing the new state, \code{x}, and whether the proposed value was accepted, logical \code{accpt}.
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' n_iter <- 100 # set to 1e3 for more complete illustration
#' draws <- matrix(0.2, nrow = n_iter, ncol = 2)
#' nAccpt <- 0L
#' pseudo <- list( list(ld = function(x) dbeta(x, 2, 2, log = TRUE),
#'                      q = function(u) qbeta(u, 2, 2)),
#'                 list(ld = function(x) dbeta(x, 2, 2, log = TRUE),
#'                      q = function(u) qbeta(u, 2, 2))
#' )
#' for (i in seq.int(2, n_iter)) {
#'  out <- imh_pseudo(draws[i - 1, ], log_target = lf, pseudo = pseudo)
#'  draws[i,] <- out$x
#'  nAccpt <- nAccpt + out$accpt
#'  cat(i, '\r')
#' }
#' nAccpt / (nrow(draws) - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
imh_pseudo <- function(x, log_target, pseudo) {
  K <- length(x)
  if (K == 1) {
    out <- imh_pseudo_univ(x = x, log_target = log_target, pseudo = pseudo, K = K)
  } else {
    out <- imh_pseudo_mv(x = x, log_target = log_target, pseudo = pseudo, K = K)
  }
  out
}

imh_pseudo_univ <- function(x, log_target, pseudo, K = K) {

  lfx0 <- log_target(x) - pseudo$ld(x)

  u1 <- runif(K)

  x1 <- pseudo$q(u1)

  lfx1 <- log_target(x1) - pseudo$ld(x1)

  lprob_accpt <- lfx1 - lfx0
  lu_accpt <- log(runif(1))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(x = x1, accpt = TRUE)
  } else {
    out <- list(x = x, accpt = FALSE)
  }
  out
}

imh_pseudo_mv <- function(x, log_target, pseudo, K = K) {

  Kpseu <- length(pseudo)

  if (Kpseu == K) {

    lfx0 <- log_target(x) - sum(sapply(1:K, function(k) pseudo[[k]]$ld(x[k])))

    u1 <- runif(K)
    x1 <- sapply(1:K, function(k) pseudo[[k]]$q(u1[k]))
    lfx1 <- log_target(x1) - sum(sapply(1:K, function(k) pseudo[[k]]$ld(x1[k])))

  } else if (Kpseu == 1) {

    lfx0 <- log_target(x) - pseudo$ld(x)
    x1 <- pseudo$r()
    lfx1 <- log_target(x1) - pseudo$ld(x1)

  } else {
    stop("imh_pseudo_mv() requires the pseudo (proposal) list to be of length 1 or
         length matching that of input vector x.")
  }

  lprob_accpt <- lfx1 - lfx0
  lu_accpt <- log(runif(1))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(x = x1, accpt = TRUE)
  } else {
    out <- list(x = x, accpt = FALSE)
  }
  out
}
