#' Independence Metropolis-Hastings
#'
#' @param x The current state (scalar or numeric vector).
#' @param log_target A function taking a scalar or numeric vector that evaluates the log-target
#'   density, returning a numeric scalar.
#' @param pseudo List specifying the pseudo-target (proposal distribution). If the list length is
#' equal to the number of dimensions in \code{x}, each element is itself a list that specifies
#' the pseudo-target independently for each corresponding dimension with functions \code{ld}
#' that evaluates the log density for that dimension,
#' and \code{q} that evaluates the quantile (inverse-CDF) function for that dimension.
#' If the dimension of \code{x} is one, then supply only the inner list
#' specifying the single pseudo-target.
#'
#' If \code{x} is a vector but a single pseudo-target is supplied, the list must
#' contain a log-density function \code{ld} that accepts a vector, and a
#' function \code{r} that takes no arguments and generates a single multivariate draw from the
#' proposal distribution.
#' @return A list containing the new state, \code{x}, and whether the proposed value was accepted, logical \code{accpt}.
#' @importFrom stats runif
#' @export
#' @examples
#' lf <- function(x) dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
#' n_iter <- 10 # set to 1e3 for more complete illustration
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
  stopifnot(
    all(is.finite(x)),
    is.function(log_target),
    is.list(pseudo)
  )

  K <- length(x)
  if (K == 1) {
    out <- imh_pseudo_univ(
      x = x,
      log_target = log_target,
      pseudo = pseudo
    )
  } else {
    out <- imh_pseudo_mv(x = x, log_target = log_target, pseudo = pseudo, K = K)
  }
  out
}

imh_pseudo_univ <- function(x, log_target, pseudo) {
  stopifnot(
    is.function(pseudo$ld),
    is.function(pseudo$q)
  )

  nEvaluations <- 0L

  lf <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_target,
      x = z,
      sampler = "imh_pseudo",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }

  ld_pseudo <- function(z, stage) {
    .eval_pseudo_logdens(
      logdens = pseudo$ld,
      x = z,
      sampler = "imh_pseudo",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = FALSE
    )
  }

  lh <- function(x, stage, allow_minus_inf = TRUE) {
    lf(x, stage = stage, allow_minus_inf = allow_minus_inf) -
      ld_pseudo(x, stage = stage)
  }

  lhx0 <- lh(x, stage = "evaluate at current state", allow_minus_inf = FALSE)

  u1 <- runif(1, min = 0.0, max = 1.0)

  x1 <- .eval_pseudo_quantile(
    q = pseudo$q,
    u = u1,
    sampler = "imh_pseudo",
    stage = "generate proposal",
    attempt = 1L,
    expected_length = 1L
  )

  lhx1 <- lh(x1, stage = "evaluate at proposal")

  lprob_accpt <- lhx1 - lhx0
  lu_accpt <- log(runif(1, min = 0.0, max = 1.0))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(x = x1, accpt = TRUE, nEvaluations = nEvaluations)
  } else {
    out <- list(x = x, accpt = FALSE, nEvaluations = nEvaluations)
  }
  out
}

imh_pseudo_mv <- function(x, log_target, pseudo, K) {
  is_joint_pseudo <- is.function(pseudo$ld) && is.function(pseudo$r)
  is_component_pseudo <-
    length(pseudo) == K &&
    all(vapply(
      pseudo,
      function(z) {
        is.list(z) && is.function(z$ld) && is.function(z$q)
      },
      logical(1)
    ))

  nEvaluations <- 0L

  lf <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_target,
      x = z,
      sampler = "imh_pseudo",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }

  if (isTRUE(is_component_pseudo)) {
    ld_pseudo <- function(z, k, stage) {
      .eval_pseudo_logdens(
        logdens = pseudo[[k]]$ld,
        x = z,
        sampler = "imh_pseudo",
        stage = stage,
        evaluation = nEvaluations,
        allow_minus_inf = FALSE
      )
    }

    lhx0 <- lf(
      x,
      stage = "evaluate at the current state",
      allow_minus_inf = FALSE
    ) -
      sum(sapply(1:K, function(k) {
        ld_pseudo(x[k], k = k, stage = "evaluate at the current state")
      }))

    u1 <- runif(K, min = 0.0, max = 1.0)

    x1 <- numeric(K)
    for (k in 1:K) {
      x1[k] <- .eval_pseudo_quantile(
        q = pseudo[[k]]$q,
        u = u1[k],
        sampler = "imh_pseudo",
        stage = sprintf("component %d proposal", k),
        attempt = 1L,
        expected_length = 1L
      )
    }

    lhx1 <- lf(
      x1,
      stage = "evaluate at proposal",
      allow_minus_inf = TRUE
    ) -
      sum(sapply(1:K, function(k) {
        ld_pseudo(x1[k], k = k, stage = "evaluate at proposal")
      }))
  } else if (isTRUE(is_joint_pseudo)) {
    ld_pseudo <- function(z, stage) {
      .eval_pseudo_logdens(
        logdens = pseudo$ld,
        x = z,
        sampler = "imh_pseudo",
        stage = stage,
        evaluation = nEvaluations,
        allow_minus_inf = FALSE
      )
    }

    lhx0 <- lf(
      x,
      stage = "evaluate at the current state",
      allow_minus_inf = FALSE
    ) -
      ld_pseudo(x, stage = "evaluate at the current state")

    x1 <- pseudo$r()

    if (length(x1) != K) {
      stop(
        "Dimension of proposal, ",
        length(x1),
        ", does not match dimension of input value, ",
        K
      )
    }

    lhx1 <- lf(
      x1,
      stage = "evaluate at proposal",
      allow_minus_inf = TRUE
    ) -
      ld_pseudo(x1, stage = "evaluate at proposal")
  } else {
    stop(
      "Pseudo-target structure passed to imh_pseudo() is not supported. 
      Multivariate sampling requires the pseudo-target (proposal) list to contain a joint log-density and generation function, 
      or have length matching that of input vector x (independent components)."
    )
  }

  lprob_accpt <- lhx1 - lhx0
  lu_accpt <- log(runif(1, min = 0.0, max = 1.0))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(x = x1, accpt = TRUE, nEvaluations = nEvaluations)
  } else {
    out <- list(x = x, accpt = FALSE, nEvaluations = nEvaluations)
  }
  out
}
