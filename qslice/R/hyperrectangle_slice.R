#' Multivariate Slice Sampler with Shrinking Hyperrectangle
#'
#' Multivariate slice sampler
#' in Algorithm 8 of Neal (2003) using the "shrinkage" procedure.
#'
#' @param x The current state (as a numeric vector).
#' @param target A function taking numeric vector that evaluates the log-target
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
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' draws <- matrix(0.2, nrow = 10e3, ncol = 2)
#' nEvaluations <- 0L
#' for (i in seq.int(2, nrow(draws))) {
#'  out <- slice_hyperrect(draws[i - 1, ], target = lf, w = c(0.5, 0.5))
#'  draws[i,] <- out$x
#'  nEvaluations <- nEvaluations + out$nEvaluations
#'  cat(i, '\r')
#' }
#' nEvaluations / (nrow(draws) - 1)
#' nEvaluations / diag(ess(draws))
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#'
slice_hyperrect <- function(x, target, w = NULL, L = NULL, R = NULL) {

  k <- length(x)

  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
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

#' Multivariate Transform Slice Sampler
#'
#' Quantile slice sampler for a random vector. The pseudo-target is specified
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
#' @importFrom stats runif
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
#' draws <- matrix(0.2, nrow = 10e3, ncol = 2)
#' draws_u <- draws
#' draws_u[1,] <- sapply(1:length(ps), function(k) ps[[k]]$p(draws[1,k]))
#' nEvaluations <- 0L
#' for (i in seq.int(2, nrow(draws))) {
#'   out <- slice_mv_transform(draws[i - 1, ], target = lf, pseudo = ps)
#'   draws[i,] <- out$x
#'   draws_u[i,] <- out$u
#'   nEvaluations <- nEvaluations + out$nEvaluations
#'   cat(i, '\r')
#' }
#' nEvaluations / (nrow(draws) - 1)
#' nEvaluations / diag(ess(draws))
#' plot(draws[,1], draws[,2], xlim = c(0, 1))
#' hist(draws[,1], freq = FALSE); curve(dbeta(x, 3, 4), col = "blue", add = TRUE)
#' hist(draws[,2], freq = FALSE); curve(dbeta(x, 5, 3), col = "blue", add = TRUE)
#' plot(draws_u[,1], draws_u[,2], xlim = c(0, 1))
#' hist(draws_u[,1], freq = FALSE)
#' hist(draws_u[,2], freq = FALSE)
#' auc(u = draws_u[,1])
#' auc(u = draws_u[,2])
slice_mv_transform <- function(x, target, pseudo) {

  K <- length(x)

  lhu <- function(u) {
    xx <- sapply(1:K, function(k) pseudo[[k]]$q(u[k]))
    target(xx) - sum(sapply(1:K, function(k) pseudo[[k]]$ld(xx[k])))
  }

  u0 <- sapply(1:K, function(k) pseudo[[k]]$p(x[k]))
  shyp <- slice_hyperrect(u0, target = lhu, w = NULL, L = NULL, R = NULL)
  x1 <- sapply(1:K, function(k) pseudo[[k]]$q(shyp$x[k]))

  list(x = x1, u = shyp$x, nEvaluations = shyp$nEvaluations)
}






#' Sequence of Conditional Pseudo-Targets from a Realization
#'
#' Given a realization of a random vector, generate a the corresponding
#' sequence of conditional pseudo-target inverse cdfs.
#' The pseudo-target is specified as
#' a sequence of growing conditional distributions.
#'
#' @param x A numeric vector of values between 0 and 1.
#' @param pseu_init A list output from \code{pseu_t_list()} describing the
#' marginal pseudo-target for \code{x[1]}.
#' @param loc_fn A function that specifies the location of a conditional
#'  pseudo-target given the elements in \code{x} that precede it.
#' @param sc_fn A function that specifies the scale of a conditional
#'  pseudo-target given the elements in \code{x} that precede it
#' @param degf A positive scalar for the degrees of freedom for all conditional
#' pseudo-targets.
#' @param lb A numeric vector (same length as \code{x}) specifying the lower
#'  bound of support for each conditional pseudo-target.
#' @param ub A numeric vector (same length as \code{x}) specifying the upper
#'  bound of support for each conditional pseudo-target.
#'
#' @return A list containing \code{x} obtained from the sequence of inverse
#' cdfs, and \code{pseudo_t_seq}, a list output from \code{pseu_t_list()}
#' describing the sequence of conditional pseudo-targets.
#'
#' @importFrom stats runif
#' @export
pseudo_t_condseq <- function(x, pseu_init, loc_fn, sc_fn, degf, lb, ub) {

  # lb and ub must have length K (even though first elements will be ignored)

  K <- length(x)
  out <- list()
  out[[1]] <- pseu_init

  for (k in 2:K) {
    out[[k]] <- pseudo_t_list(loc = loc_fn(x[1:(k-1)]), sc = sc_fn(x[1:(k-1)]),
                              degf = degf, lb = lb[k], ub = ub[k])
  }

  out
}


#' Inverse Transform from Sequence of Conditional Pseudo-Targets
#'
#' Given a vector of from a unit hypercube, map to the original vector using
#' a sequence of conditional pseudo-target inverse cdfs.
#' The pseudo-target is specified as
#' a sequence of growing conditional distributions.
#'
#' @param u A numeric vector of values between 0 and 1.
#' @param pseu_init A list output from \code{pseu_t_list()} describing the
#' marginal pseudo-target for \code{x[1]}.
#' @param loc_fn A function that specifies the location of a conditional
#'  pseudo-target given the elements in \code{x} that precede it.
#' @param sc_fn A function that specifies the scale of a conditional
#'  pseudo-target given the elements in \code{x} that precede it
#' @param degf A positive scalar for the degrees of freedom for all conditional
#' pseudo-targets.
#' @param lb A numeric vector (same length as \code{x}) specifying the lower
#'  bound of support for each conditional pseudo-target.
#' @param ub A numeric vector (same length as \code{x}) specifying the upper
#'  bound of support for each conditional pseudo-target.
#'
#' @return A list containing \code{x} obtained from the sequence of inverse
#' cdfs, and \code{pseudo_t_seq}, a list output from \code{pseu_t_list()}
#' describing the sequence of conditional pseudo-targets.
#'
#' @importFrom stats runif
#' @export
pseudo_t_condseq_XfromU <- function(u, pseu_init, loc_fn, sc_fn, degf, lb, ub) {

  # lb and ub must have length K (even though first elements will be ignored)

  K <- length(u)
  out <- list()
  out[[1]] <- pseu_init

  x <- numeric(K)

  x[1] <- pseu_init$q(u[1])

  for (k in 2:K) {
    out[[k]] <- pseudo_t_list(loc = loc_fn(x[1:(k-1)]), sc = sc_fn(x[1:(k-1)]),
                              degf = degf, lb = lb[k], ub = ub[k])
    x[k] <- out[[k]]$q(u[k])
  }

  list(x = x, pseudo_t_seq = out)
}




#' Multivariate Transform Slice Sampler from a Sequence of Conditional Pseudo-Targets
#'
#' Quantile slice sampler for a random vector. The pseudo-target is specified as
#' a sequence of growing conditional distributions.
#'
#' @inherit slice_hyperrect
#' @param pseudo_control A list with \code{pseu_init}, a list output from
#' \code{pseu_t_list()} describing the marginal pseudo-target for \code{x[1]};
#' \code{loc_fn}, a function that specifies the location of a conditional
#'  pseudo-target given the elements in \code{x} that precede it; \code{sc_fn},
#'  a function that specifies the scale of a conditional
#'  pseudo-target given the elements in \code{x} that precede it; \code{degf}, a
#'  positive scalar for the degrees of freedom for all conditional pseudo-targets;
#'  \code{lb}, a numeric vector (same length as \code{x}) specifying the lower
#'  bound of support for each conditional pseudo-target;
#'  \code{ub}, a numeric vector (same length as \code{x}) specifying the upper
#'  bound of support for each conditional pseudo-target.
#'
#' @return A list containing three elements: "x" is the new state, "u" is a vector
#'   of values of the sequence of conditional CDFs of the psuedo-targets associated
#'   with the returned value, and "nEvaluations is the number of evaluations of the
#'   target function used to obtain the new state.
#'
#' @importFrom stats runif
#' @seealso \pkg{\link{coda}}
#' @export
#' @examples
#' # Funnel distribution from Neal (2003).
#' K <- 10
#' n_iter <- 10e3 # MCMC iterations
#' n <- 10e3 # number of iid samples from the target
#' Y <- matrix(NA, nrow = n, ncol = K) # iid samples from the target
#' Y[,1] <- rnorm(n, 0.0, 3.0)
#' for (i in 1:n) {
#'   Y[i, 2:K] <- rnorm(K-1, 0.0, exp(0.5*Y[i,1]))
#' }
#'
#' ltarget <- function(x) {
#' dnorm(x[1], 0.0, 3.0, log = TRUE) +
#'   sum(dnorm(x[2:K], 0.0, exp(0.5*x[1]), log = TRUE))
#' }
#'
#' pseudo_control <- list(
#'
#'   loc_fn = function(x) {
#'     0.0
#'   },
#'
#'   sc_fn = function(x) {
#'     if (is.null(x)) {
#'       out <- 3.0
#'     } else {
#'       out <- exp(0.5*x[1])
#'     }
#'     out
#'   },
#'
#'   pseu_init = pseudo_t_list(loc = 0.0, sc = 3.0, degf = 5.0,
#'                             lb = -Inf, ub = Inf),
#'
#'   degf = 5.0,
#'   lb = rep(-Inf, K),
#'   ub = rep(Inf, K)
#' )
#'
#' x0 <- runif(K)
#' draws <- matrix(rep(x0, n_iter + 1), nrow = n_iter + 1, byrow = TRUE)
#' draws_u <- matrix(rep(x0, n_iter), nrow = n_iter, byrow = TRUE)
#' n_eval <- 0
#' for (i in 2:(n_iter + 1)) {
#'   tmp <- slice_mv_transform_seq(draws[i-1,],
#'                                 target = ltarget,
#'                                 pseudo_control = pseudo_control)
#'   draws[i,] <- tmp$x
#'   draws_u[i-1,] <- tmp$u
#'   n_eval <- n_eval + tmp$nEvaluations
#' }
#'
#' if (requireNamespace("coda", quietly = TRUE)) {
#'   (es <- coda::effectiveSize(coda::as.mcmc(draws)))
#'   mean(es)
#' }
#'
#' n_eval / n_iter
#' sapply(1:K, function (k) auc(u = draws_u[,k]))
#'
#' hist(draws_u[,1])
#' plot(draws[,1], draws[,2])
#' points(Y[,1], Y[,2], col = "blue", cex = 0.5)
#'
slice_mv_transform_seq <- function(x, target, pseudo_control) {

  # target is log target only without pseudo
  # pseudo_control is a list with: pseu_init, loc_fn, sc_fn, degf, lb, ub

  K <- length(x)

  # Step 1
  tmp_seq <- pseudo_t_condseq(x = x, pseu_init = pseudo_control$pseu_init,
                              loc_fn = pseudo_control$loc_fn,
                              sc_fn = pseudo_control$sc_fn,
                              degf = pseudo_control$degf,
                              lb = pseudo_control$lb,
                              ub = pseudo_control$ub)

  fx <- target(x) - sum(sapply(1:K, function(k) tmp_seq[[k]]$ld(x[k])))
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
    tmp <- pseudo_t_condseq_XfromU(u = u1, pseu_init = pseudo_control$pseu_init,
                                   loc_fn = pseudo_control$loc_fn,
                                   sc_fn = pseudo_control$sc_fn,
                                   degf = pseudo_control$degf,
                                   lb = pseudo_control$lb,
                                   ub = pseudo_control$ub)
    x1 <- tmp$x
    tmp_seq <- tmp$pseudo_t_seq

    fx1 <- target(x1) - sum(sapply(1:K, function(k) tmp_seq[[k]]$ld(x1[k])))
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
