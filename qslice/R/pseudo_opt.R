#' Optimal pseudo-target for a given target
#'
#' Find an optimal pseudo-target in a specified family to approximate
#' the given (unnormalized) target (Heiner et al., 2024+). Optimize over the selected utility function.
#'
#' Optionally supply samples from the target distribution.
#'
#' @inherit utility_pseudo
#' @param samples Optional numeric vector providing samples from the target distribution
#' (for use as alternative to \code{log_target}).
#' @param type String specifying the input type. One of "function", "samples", or "grid".
#' Default is to use "samples".
#'
#' Use of "function" requires specification of \code{log_target}.
#'
#' Use of "samples" requires specification of \code{samples}.
#' @param family String specifying the family of distributions for the pseudo-target.
#' Can be any of the families accepted by \link[qslice]{pseudo_list}.
#' @param degf Numeric vector of degrees of freedom values to try (only if \code{family = "t"}.
#' Defaults to \code{c(1, 5, 20)}.
#' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
#' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
#' @param nbins Positive integer specifying the number of histogram bins if using "samples" or "grid".
#' Defaults to 100.
#' @param tol_opt Positive numeric scalar that passes to \code{reltol} in the call
#' to \link[stats]{optim}. Defaults to \code{1.0e-6}.
#' @param tol_int Positive numeric scalar that passes to \code{abs.tol} in the call to \link[stats]{integrate}.
#' Defaults to \code{1.0e-3}.
#' @param verbose Logical for whether to print intermediate steps of optimization.
#' Defaults to \code{FALSE}.
#' @returns A list with named components:
#'
#'  \code{pseudo}: a list with functions corresponding to the selected pseudo-target;
#'  output of \link[qslice]{pseudo_list}.
#'
#'  \code{utility}: value of the utility function using the selected pseudo-target;
#'  output of \link[qslice]{utility_pseudo}.
#'
#'  \code{utility_type}: repeats the input specifying the utility type.
#'
#'  \code{opt}: output of \link[stats]{optim}.
#'
#'  Other outputs repeating inputs.
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*
#'
#' @importFrom stats sd optim
#' @importFrom stats sd
#'
#' @export
#' @examples
#' (pseu <- pseudo_opt(samples = rnorm(1e3), type = "samples",
#'                family = "t", utility_type = "AUC",
#'                nbins = 10, plot = TRUE,
#'                verbose = FALSE))
#' oldpar <- par(mfrow = c(1,2))
#' (pseu <- pseudo_opt(log_target = function(x) dnorm(x, log = TRUE),
#'                 type = "function",
#'                 family = "logistic", utility_type = "AUC",
#'                 nbins = 100, plot = TRUE,
#'                 verbose = FALSE))
#' (pseu <- pseudo_opt(log_target = function(x) dbeta(x, 4, 2, log = TRUE),
#'                 lb = 0, ub = 1,
#'                 type = "function",
#'                 family = "cauchy", utility_type = "AUC",
#'                 nbins = 30, plot = TRUE,
#'                 verbose = FALSE))
#' par(oldpar)
#'
pseudo_opt <- function(log_target = NULL,
                       samples = NULL,
                       type = "samples", # one of "samples", "function"
                       family = "t",
                       degf = c(1, 5, 20),
                       lb = -Inf, ub = Inf,
                       utility_type = "AUC",
                       nbins = 100,
                       tol_opt = 1.0e-6, tol_int = 1.0e-3,
                       plot = TRUE,
                       verbose = FALSE) {

  if (family == "beta") stop("Optimazation for beta pseudo-targets not implemented.")

  x <- seq(1.0e-6, 1.0 - 1.0e-6, length = nbins) # trouble if it doesn't reach far enough into tails

  if (is.null(samples)) {
    inits <- c(loc = 0.5, sc = 2.0)
  } else {
    inits <- c(loc = mean(samples), sc = sd(samples))
  }

  get_util <- function(pars, type, x, log_target, samples, family,
                       degf, lb, ub, utility_type, nbins, verbose) {

    loc <- pars[1]
    sc <- pars[2]

    if (sc <= 0.0) {

      out <- -1.0

    } else {

      pseu <- pseudo_list(family = family,
                          params = list(loc = loc, sc = sc, degf = degf),
                          lb = lb, ub = ub)

      if (verbose) cat("trying", pseu$t, "\n")

      # this function is found in utility_shrinkslice.R; output is the same as utility_shrinkslice()
      out <- utility_pseudo(pseudo = pseu,
                            log_target = log_target,
                            samples = samples,
                            type = type, x = x,
                            nbins = nbins,
                            plot = FALSE,
                            utility_type = utility_type,
                            tol_int = tol_int)

      if (verbose) {
        print(out)
        cat("\n")
      }

      if (out > 1.0) { # then the result isn't viable. Penalize.
        out <- -1.0
      }

    }

    out
  }

  opt <- list()

  if (family != "t") {
    degf = NA
  }

  for (i in 1:length(degf)) {
    opt[[i]] <- optim(inits, get_util, control = list(fnscale=-1, reltol = tol_opt),
                      type = type,
                      x = x, log_target = log_target, samples = samples,
                      family = family,
                      degf = degf[i],
                      lb = lb, ub = ub,
                      utility_type = utility_type,
                      nbins = nbins,
                      verbose = verbose)
   }

  use_indx <- which.max(sapply(opt, function(obj) obj$value))

  pseu <- pseudo_list(family = family,
                      params = list(loc = unname(opt[[use_indx]]$par[1]),
                                    sc = unname(opt[[use_indx]]$par[2]),
                                    degf = degf[use_indx]),
                      lb = lb, ub = ub)

  util <- utility_pseudo(pseudo = pseu,
                         log_target = log_target,
                         samples = samples,
                         type = type, x = x,
                         nbins = nbins,
                         plot = plot,
                         utility_type = utility_type,
                         tol_int = tol_int)

  list(pseudo = pseu, utility = util, utility_type = utility_type,
       opt = opt[[use_indx]],
       nbins = nbins, tol_int = tol_int, tol_opt = tol_opt)
}
