#' Multivariate Slice Sampler with Shrinking Hyperrectangle
#'
#' Multivariate slice sampler
#' in Algorithm 8 of Neal (2003) using the "shrinkage" procedure.
#'
#' @param x The current state (as a numeric vector).
#' @param log_target A function taking numeric vector that evaluates the log-target
#'   density, returning a numeric scalar.
#' @param w A numeric vector tuning the algorithm which gives the typical slice
#'   width in each dimension. This is a main tuning parameter of the algorithm.
#'   If \code{NULL}, the sampler begins shrinking from the supplied boundaries (should,
#'   correspond with the support).
#' @param L Numeric vector giving the lower boundary of support in each dimension.
#' @param R Numeric vector giving the upper boundary of support in each dimension.
#' Will be used if \code{w} is null. If all of \code{L}, \code{R}, and \code{w}
#' are null, then the boundaries default to those of the unit hypercube.
#'
#' @return A list contains two elements: "x" is the new state and "nEvaluations"
#'   is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#' @references
#' Neal, R. M. (2003), "Slice sampling," *The Annals of Statistics*, 31, 705-767. \doi{https://doi.org/10.1214/aos/1056562461}
#'
#' @importFrom stats runif
#' @importFrom graphics curve points segments text
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' n_iter <- 10 # set to 1e4 for more complete illustration
#' draws <- matrix(0.2, nrow = n_iter, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, n_iter)) {
#'  out <- slice_hyperrect(draws[i - 1, ], log_target = lf, w = c(0.5, 0.5))
#'  draws[i,] <- out$x
#'  nEvaluations <- nEvaluations + out$nEvaluations
#'  cat(i, '\r')
#' }
#' nEvaluations / (nrow(draws) - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#'
slice_hyperrect <- function(x, log_target, w = NULL, L = NULL, R = NULL) {

  k <- length(x)

  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x)
  }

  # Step 1
  fx <- f(x)
  stopifnot(fx > -Inf)
  y <- log(runif(1)) + fx

  # Step 2 (Initialize hyperrectangle)

  if (!all(is.null(w))) {
    L <- x - runif(k) * w
    R <- L + w
  } else if (all(is.null(c(L, R)))) {
    L <- rep(0.0, k)
    R <- rep(1.0, k)
  }

  # Step 3 ("Shrinkage" procedure)
  repeat {

    x1 <- L + runif(k) * (R - L)

    if (y < f(x1)) {
      return(list(x = x1, nEvaluations = nEvaluations))
    }

    shrinkL <- (x1 < x)

    if (any(shrinkL)) {
      indx <- which(shrinkL)
      L[indx] <- x1[indx]
    }

    if (any(!shrinkL)) {
      indx <- which(!shrinkL)
      R[indx] <- x1[indx]
    }

  }
}

#' Multivariate Quantile Slice Sampler
#'
#' Quantile slice sampler for a random vector (Heiner et al., 2024+). The pseudo-target is specified
#' through independent univariate distributions.
#'
#' @inherit slice_hyperrect
#' @param pseudo List of length equal to the number of dimensions in \code{x}. Each element is itself a list that specifies
#' the pseudo-target for the corresponding dimension with functions \code{ld} that evaluates the log density,
#' \code{p} that evaluates the CDF, and \code{q} that evaluates the quantile (inverse-CDF) function.
#' @return A list containing three elements: "x" is the new state, "u" is the
#'   value of the CDF of the psuedo-target associated with the returned value,
#'   inverse CDF method, and "nEvaluations is the number of evaluations of the
#'   target function used to obtain the new state.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*
#'
#' @importFrom stats runif
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' ps_shsc <- list(c(2, 2), c(2, 1))
#' ps <- list(
#'   list(ld = function(x) dbeta(x, ps_shsc[[1]][1], ps_shsc[[1]][2], log = TRUE),
#'        p = function(x) pbeta(x, ps_shsc[[1]][1], ps_shsc[[1]][2]),
#'        q = function(x) qbeta(x, ps_shsc[[1]][1], ps_shsc[[1]][2]) ),
#'   list(ld = function(x) dbeta(x, ps_shsc[[2]][1], ps_shsc[[2]][2], log = TRUE),
#'        p = function(x) pbeta(x, ps_shsc[[2]][1], ps_shsc[[2]][2]),
#'        q = function(x) qbeta(x, ps_shsc[[2]][1], ps_shsc[[2]][2]) )
#'   )
#' n_iter <- 10 # set to 1e4 for more complete illustration
#' draws <- matrix(0.2, nrow = n_iter, ncol = 2)
#' draws_u <- draws
#' draws_u[1,] <- sapply(1:length(ps), function(k) ps[[k]]$p(draws[1,k]))
#' nEvaluations <- 0L
#' for (i in seq.int(2, n_iter)) {
#'   out <- slice_quantile_mv(draws[i - 1, ], log_target = lf, pseudo = ps)
#'   draws[i,] <- out$x
#'   draws_u[i,] <- out$u
#'   nEvaluations <- nEvaluations + out$nEvaluations
#'   cat(i, '\r')
#' }
#' nEvaluations / (nrow(draws) - 1)
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#' plot(draws_u[,1], draws_u[,2], xlim = c(0, 1))
#' hist(draws_u[,1], freq = FALSE)
#' hist(draws_u[,2], freq = FALSE)
#' auc(u = draws_u[,1])
#' auc(u = draws_u[,2])
slice_quantile_mv <- function(x, log_target, pseudo) {

  K <- length(x)

  lhu <- function(u) {
    xx <- sapply(1:K, function(k) pseudo[[k]]$q(u[k]))
    log_target(xx) - sum(sapply(1:K, function(k) pseudo[[k]]$ld(xx[k])))
  }

  u0 <- sapply(1:K, function(k) pseudo[[k]]$p(x[k]))
  shyp <- slice_hyperrect(u0, log_target = lhu, w = NULL, L = NULL, R = NULL)
  x1 <- sapply(1:K, function(k) pseudo[[k]]$q(shyp$x[k]))

  list(x = x1, u = shyp$x, nEvaluations = shyp$nEvaluations)
}






#' Sequence of conditional pseudo-targets from a realization
#'
#' Given a realization of a random vector, generate a the corresponding
#' sequence of conditional pseudo-target inverse CDFs (Heiner et al., 2024+).
#' The pseudo-target is specified as
#' a sequence of growing conditional distributions.
#'
#' See the documentation for \link[qslice]{slice_quantile_mv_seq} for examples.
#'
#' @param x A numeric vector of values between 0 and 1.
#' @param pseudo_init A list output from \link[qslice]{pseudo_list} describing the
#' marginal pseudo-target for \code{x[1]}. All subsequent pseudo-targets will
#' resemble \code{pseudo_init} with exception of different location and scale parameters.
#' @param loc_fn A function that specifies the location of a conditional
#'  pseudo-target given the elements in \code{x} that precede it.
#' @param sc_fn A function that specifies the scale of a conditional
#'  pseudo-target given the elements in \code{x} that precede it
#' @param lb A numeric vector (same length as \code{x}) specifying the lower
#'  bound of support for each conditional pseudo-target.
#' @param ub A numeric vector (same length as \code{x}) specifying the upper
#'  bound of support for each conditional pseudo-target.
#'
#' @return A list containing a sequence of pseudo-targets, each from \link[qslice]{pseudo_list}.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*
#'
#' @importFrom stats runif
#'
#' @export
#' @example man/examples/pseudo_sequential.R
#'
pseudo_condseq <- function(x, pseudo_init, loc_fn, sc_fn, lb, ub) {

  # lb and ub must have length K (even though first elements will be ignored)

  K <- length(x)
  out <- list()
  out[[1]] <- pseudo_init

  family <- pseudo_init$family
  params_now <- pseudo_init$params

  for (k in 2:K) {

    params_now$loc <- loc_fn(x[1:(k-1)])
    params_now$sc <- sc_fn(x[1:(k-1)])

    out[[k]] <- pseudo_list(family = family,
                            params = params_now,
                            lb = lb[k], ub = ub[k])
  }

  out
}


#' Inverse transform from sequence of conditional pseudo-targets
#'
#' Given a vector of from a unit hypercube, map to the original (back-transformed)
#' vector using a sequence of conditional pseudo-target inverse CDFs.
#' The pseudo-target is specified as
#' a sequence of growing conditional distributions.
#'
#' See the documentation for \link[qslice]{slice_quantile_mv_seq} for examples.
#'
#' @param u A numeric vector of values between 0 and 1.
#' @param pseudo_init A list output from \link[qslice]{pseudo_list} describing the
#' marginal pseudo-target for \code{x[1]}.
#' @param loc_fn A function that specifies the location of a conditional
#'  pseudo-target given the elements in \code{x} that precede it.
#' @param sc_fn A function that specifies the scale of a conditional
#'  pseudo-target given the elements in \code{x} that precede it
#' @param lb A numeric vector (same length as \code{x}) specifying the lower
#'  bound of support for each conditional pseudo-target.
#' @param ub A numeric vector (same length as \code{x}) specifying the upper
#'  bound of support for each conditional pseudo-target.
#'
#' @return A list containing \code{x} obtained from the sequence of inverse
#' CDFs, and \code{pseudo_seq}, a list of the corresponding sequential
#' pseudo-targets output from \link[qslice]{pseudo_list}.
#'
#' @export
#' @example man/examples/pseudo_sequential.R
#'
pseudo_condseq_XfromU <- function(u, pseudo_init, loc_fn, sc_fn, lb, ub) {

  # lb and ub must have length K (even though first elements will be ignored)

  K <- length(u)
  out <- list()
  out[[1]] <- pseudo_init

  x <- numeric(K)
  x[1] <- pseudo_init$q(u[1])
  family <- pseudo_init$family
  params_now <- pseudo_init$params

  for (k in 2:K) {

    params_now$loc <- loc_fn(x[1:(k-1)])
    params_now$sc <- sc_fn(x[1:(k-1)])

    out[[k]] <- pseudo_list(family = family,
                            params = params_now,
                            lb = lb[k], ub = ub[k])
    x[k] <- out[[k]]$q(u[k])
  }

  list(x = x, pseudo_seq = out)
}



#' Multivariate Quantile Slice Sampler from a sequence of conditional pseudo-targets
#'
#' Quantile slice sampler for a random vector (Heiner et al., 2024+). The pseudo-target is specified as
#' a sequence of growing conditional distributions.
#'
#' @inherit slice_hyperrect
#' @param pseudo_control A list with
#'
#' \code{pseudo_init}, a list output from
#' \link[qslice]{pseudo_list} describing the marginal pseudo-target for \code{x[1]}.
#' Attributes of \code{pseudo_init} will be used in subsequent pseudo-targets,
#' except for location and scale parameters.
#'
#' \code{loc_fn}, a function that specifies the location of a conditional
#'  pseudo-target given the elements in \code{x} that precede it.
#'
#'  \code{sc_fn}, a function that specifies the scale of a conditional
#'  pseudo-target given the elements in \code{x} that precede it.
#'
#'  \code{lb}, a numeric vector (same length as \code{x}) specifying the lower
#'  bound of support for each conditional pseudo-target.
#'
#'  \code{ub}, a numeric vector (same length as \code{x}) specifying the upper
#'  bound of support for each conditional pseudo-target.
#'
#' @return A list containing three elements: "x" is the new state, "u" is a vector
#'   of values of the sequence of conditional CDFs of the psuedo-targets associated
#'   with the returned value, and "nEvaluations is the number of evaluations of the
#'   target function used to obtain the new state.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*
#'
#' @importFrom stats runif
#'
#' @export
#' @example man/examples/pseudo_sequential.R
#'
slice_quantile_mv_seq <- function(x, log_target, pseudo_control) {

  # target is log target only without pseudo

  K <- length(x)

  # Step 1
  tmp_seq <- pseudo_condseq(x = x,
                            pseudo_init = pseudo_control$pseudo_init,
                            loc_fn = pseudo_control$loc_fn,
                            sc_fn = pseudo_control$sc_fn,
                            lb = pseudo_control$lb,
                            ub = pseudo_control$ub)

  fx <- log_target(x) - sum(sapply(1:K, function(k) tmp_seq[[k]]$ld(x[k])))
  nEvaluations <- 1
  stopifnot(fx > -Inf)

  y <- log(runif(1)) + fx

  # Step 2 (Initialize hypercube)

  u0 <- sapply(1:K, function(k) tmp_seq[[k]]$p(x[k]))

  L <- rep(0.0, K)
  R <- rep(1.0, K)

  # Step 3 ("Shrinkage" procedure)

  repeat {

    u1 <- L + runif(K) * (R - L)
    tmp <- pseudo_condseq_XfromU(u = u1,
                                 pseudo_init = pseudo_control$pseudo_init,
                                 loc_fn = pseudo_control$loc_fn,
                                 sc_fn = pseudo_control$sc_fn,
                                 lb = pseudo_control$lb,
                                 ub = pseudo_control$ub)
    x1 <- tmp$x
    tmp_seq <- tmp$pseudo_seq

    fx1 <- log_target(x1) - sum(sapply(1:K, function(k) tmp_seq[[k]]$ld(x1[k])))
    nEvaluations <- nEvaluations + 1

    if (y < fx1) {
      return(list(x = x1, u = u1, nEvaluations = nEvaluations))
    }

    shrinkL <- (u1 < u0)

    if (any(shrinkL)) {
      indx <- which(shrinkL)
      L[indx] <- u1[indx]
    }

    if (any(!shrinkL)) {
      indx <- which(!shrinkL)
      R[indx] <- u1[indx]
    }

  }
}
