#' @export
imh_pseudo_univ <- function(x, ltarget, pseudo) {

  lfx0 <- ltarget(x) - pseu$ld(x)

  u1 <- runif(K)

  x1 <- pseu$q(u1)

  lfx1 <- ltarget(x1) - pseu$ld(x1)

  lprob_accpt <- lfx1 - lfx0
  lu_accpt <- log(runif(1))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(x = x1, accpt = TRUE)
  } else {
    out <- list(x = x, accpt = FALSE)
  }
  out
}

#' @export
imh_pseudo_mv <- function(x, ltarget, pseudo) {

  K <- length(x)

  lfx0 <- ltarget(x) - sum(sapply(1:K, function(k) pseudo[[k]]$ld(x[k])))

  u1 <- runif(K)

  x1 <- sapply(1:K, function(k) pseudo[[k]]$q(u1[k]))

  lfx1 <- ltarget(x1) - sum(sapply(1:K, function(k) pseudo[[k]]$ld(x1[k])))

  lprob_accpt <- lfx1 - lfx0
  lu_accpt <- log(runif(1))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(x = x1, accpt = TRUE)
  } else {
    out <- list(x = x, accpt = FALSE)
  }
  out
}

#' Independence Metropolis-Hastings Sampler
#'
#' @param x The current state (scalar or numeric vector).
#' @param ltarget A function taking a scalar or numeric vector that evaluates the log-target
#'   density, returning a numeric scalar.
#' @param pseudo List of length equal to the number of dimensions in \code{x}. Each element is itself a list that specifies
#' the pseudo-target for the corresponding dimension with functions \code{ld} that evaluates the log density
#' and \code{q} that evaluates the quantile (inverse-CDF) function. If the dimension of \code{x} is
#' one, then supply only the inner list specifying the single pseudo-target.
#' @return A numeric vector with the new state.
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' draws <- matrix(0.2, nrow = 10e3, ncol = 2)
#' nAccpt <- 0L
#' pseudo <- list( list(ld = function(x) dbeta(x, 2, 2, log = TRUE),
#'                      q = function(u) qbeta(u, 2, 2)),
#'                 list(ld = function(x) dbeta(x, 2, 2, log = TRUE),
#'                      q = function(u) qbeta(u, 2, 2))
#' )
#' for (i in seq.int(2, nrow(draws))) {
#'  out <- imh_pseudo(draws[i - 1, ], ltarget = lf, pseudo = pseudo)
#'  draws[i,] <- out$x
#'  nAccpt <- nAccpt + out$accpt
#'  cat(i, '\r')
#' }
#' nAccpt / (nrow(draws) - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
imh_pseudo <- function(x, ltarget, pseudo) {
  K <- length(x)
  if (K == 1) {
    out <- imh_pseudo_univ(x = x, ltarget = ltarget, pseudo = pseudo)
  } else {
    out <- imh_pseudo_mv(x = x, ltarget = ltarget, pseudo = pseudo)
  }
  out
}