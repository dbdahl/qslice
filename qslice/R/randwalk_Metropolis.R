.rwm_norm_univ_eval <- function(x, log_target, cand_sd, support) {
  stopifnot(
    length(x) == 1L,
    is.finite(x),
    is.function(log_target),
    length(cand_sd) == 1L,
    is.finite(cand_sd),
    cand_sd > 0.0,
    length(support) == 2L,
    all(!is.na(support)),
    support[2] > support[1]
  )

  nEvaluations <- 0L

  lf <- function(z, stage, allow_minus_inf = TRUE) {
    nEvaluations <<- nEvaluations + 1L

    .eval_log_target(
      log_target = log_target,
      x = z,
      sampler = "rwm_norm_univ",
      stage = stage,
      evaluation = nEvaluations,
      allow_minus_inf = allow_minus_inf
    )
  }

  x_cand <- rnorm(1, mean = x, sd = cand_sd)

  if (x_cand >= support[1] && x_cand <= support[2]) {
    lprob_accpt <- lf(
      x_cand,
      stage = "evaluate at candidate",
      allow_minus_inf = TRUE
    ) -
      lf(x, stage = "evaluate at old value", allow_minus_inf = FALSE)

    lu <- log(runif(1, min = 0.0, max = 1.0))
    accept <- isTRUE(lu < lprob_accpt)
  } else {
    accept <- FALSE
  }

  x_out <- ifelse(accept, x_cand, x)

  list(x = x_out, accept = accept, nEvaluations = nEvaluations)
}
.rwm_norm_univ_eval_fast <- function(x, log_target, cand_sd, support) {
  nEvaluations <- 0L
  lf <- function(z) {
    nEvaluations <<- nEvaluations + 1L
    log_target(z)
  }

  x_cand <- rnorm(1, mean = x, sd = cand_sd)

  if (x_cand >= support[1] && x_cand <= support[2]) {
    lprob_accpt <- lf(x_cand) - lf(x)
    lu <- log(runif(1, min = 0.0, max = 1.0))
    accept <- isTRUE(lu < lprob_accpt)
  } else {
    accept <- FALSE
  }

  x_out <- ifelse(accept, x_cand, x)

  list(x = x_out, accept = accept, nEvaluations = nEvaluations)
}
#' Univariate Random Walk Metropolis with Normal Proposal
#'
#' Single update for the univariate random-walk Metropolis sampler with a normal proposal distribution.
#'
#' @param x The current state (as a numeric scalar).
#' @param log_target A function taking numeric scalar that evaluates the
#' (potentially unnormalized) log-target density, returning a numeric scalar.
#' @param cand_sd A positive numeric scalar tuning the algorithm. Gives the standard deviation of the normal random-walk proposal.
#' @param support Numeric vector of length two giving the lower and upper boundaries of the support for the target distribution.
#' Default is \code{c(-Inf, Inf)}, in which case the target density will always be evaluated at the proposed value.
#' @param diagnostics Logical, indicating whether the sampler should be run in diagnostic mode for
#' informative error handling. When \code{diagnostics} is \code{NULL}, it is reset to match the
#' value of \code{getOption("qslice.diagnostic")} which defaults to \code{FALSE} if the option is not
#' set. Default for \code{diagnostics} is \code{NULL}.
#'
#' @return A list with three elements:
#'
#' \code{x} is the new state.
#'
#' \code{accept} Logical indicating whether the proposal was accepted.
#'
#' \code{nEvaluations} is the number of evaluations of the target function used to obtain the new
#'   state.
#'
#'
#' @importFrom stats runif rnorm
#'
#' @export
#' @examples
#' lf <- function(x) dbeta(x, 3, 4, log = TRUE)
#' draws <- numeric(10) + 0.5 # set to numeric(1e3) for more complete illustration
#' nEvaluations <- 0L
#' nAccept <- 0L
#' for (i in seq.int(2, length(draws))) {
#'   out <- rwm_norm(draws[i - 1], log_target = lf, cand_sd = 0.3, support = c(0.0, 1.0))
#'   draws[i] <- out$x
#'   nEvaluations <- nEvaluations + out$nEvaluations
#'   nAccept <- nAccept + out$accept
#' }
#' nEvaluations / (length(draws) - 1)
#' nAccept / (length(draws) - 1)
#' plot(density(draws), xlim = c(0, 1))
#' curve(exp(lf(x)), 0, 1, col = "blue", add = TRUE)
#'
rwm_norm <- function(
  x,
  log_target,
  cand_sd,
  support = c(-Inf, Inf),
  diagnostics = NULL
) {
  diagnostics <- .resolve_diagnostics(diagnostics)

  if (isTRUE(diagnostics)) {
    .rwm_norm_univ_eval(
      x = x,
      log_target = log_target,
      cand_sd = cand_sd,
      support = support
    )
  } else {
    tryCatch(
      .rwm_norm_univ_eval_fast(
        x = x,
        log_target = log_target,
        cand_sd = cand_sd,
        support = support
      ),
      error = function(e) {
        message(
          "\nFor detailed diagnostics, rerun after setting:\n",
          "  options(qslice.diagnostics = TRUE)\n",
          "or call the sampler with argument diagnostics = TRUE."
        )

        stop(e)
      }
    )
  }
}
