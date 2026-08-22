## degf deprecation bridge
.resolve_df <- function(
  df,
  degf,
  default = NULL,
  caller
) {
  if (!is.null(df) && !is.null(degf)) {
    stop(
      caller,
      "(): supply only `df`; `degf` is deprecated.",
      call. = FALSE
    )
  }

  if (!is.null(degf)) {
    warning(
      caller,
      "(): `degf` is deprecated; use `df` instead.",
      call. = FALSE
    )
    return(degf)
  }

  if (is.null(df)) {
    if (is.null(default)) {
      stop(
        caller,
        "(): `df` must be supplied.",
        call. = FALSE
      )
    }

    return(default)
  }

  df
}

#' Optimal pseudo-target for a given target
#'
#' Find an optimal pseudo-target in a specified family to approximate
#' the given (unnormalized) target (Heiner et al., 2026+). Optimize over the selected utility function.
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
#' @param df Numeric vector of degrees of freedom values to try (only if \code{family = "t"}).
#' Defaults to \code{c(1, 5, 20)}.
#' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
#' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
#' @param nbins Positive integer specifying the number of histogram bins if using "samples" or "grid".
#' Defaults to 100.
#' @param tol_opt Positive numeric scalar that passes to \code{reltol} in the call
#' to \link[stats]{optim}. Defaults to \code{1.0e-6}.
#' @param utility_options List containing the following. \code{tol_int}: Positive numeric scalar that passes to \code{abs.tol} in the call to \link[stats]{integrate} used for type = "function".
#' \code{n_opt}: Positive integer giving the number of equally spaced initialization points for numerical optimization (\link[stats]{optim}) in AUC with type = "function".
#' \code{interval_opt}: Numerical vector of length two giving the interval for numerical optimization (\link[stats]{optim}) in AUC with type = "function".
#'
#' @param verbose Logical for whether to print intermediate steps of optimization.
#' Defaults to \code{FALSE}.
#' @param degf (deprecated) Numeric vector of degrees of freedom values to try (only if \code{family = "t"}). Use \code{df} instead.
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
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2026+), "Quantile Slice Sampling," *arXiv preprint arXiv:2407.12608* \doi{https://doi.org/10.48550/arXiv.2407.12608}
#' @importFrom stats sd optim
#' @importFrom stats sd
#'
#' @export
#' @examples
#' oldpar <- par(mfrow = c(1,2))
#' (pseu <- pseudo_opt(samples = rnorm(1e3), type = "samples",
#'                family = "t", utility_type = "AUC",
#'                nbins = 10, plot = TRUE,
#'                verbose = FALSE))
#' (pseu <- pseudo_opt(log_target = function(x) dnorm(x, log = TRUE),
#'                 type = "grid",
#'                 family = "logistic", utility_type = "MSW",
#'                 nbins = 30, plot = TRUE,
#'                 verbose = FALSE))
#' (pseu <- pseudo_opt(log_target = function(x) dbeta(x, 4, 2, log = TRUE),
#'                 lb = 0, ub = 1,
#'                 type = "function",
#'                 family = "cauchy", utility_type = "AUC",
#'                 nbins = 30, plot = TRUE,
#'                 verbose = FALSE))
#' par(oldpar)
#'
pseudo_opt <- function(
  log_target = NULL,
  samples = NULL,
  type = "samples", # one of "function", "samples", or "grid"
  family = "t",
  df = NULL,
  lb = -Inf,
  ub = Inf,
  utility_type = "AUC",
  nbins = 100,
  tol_opt = 1.0e-6,
  utility_options = list(
    tol_int = 1.0e-3, # integration accuracy
    n_opt = 10, # number of inits for numerical optimization in AUC with type = "integration"
    interval_opt = c(0.000001, 0.999999) # interval for numerical optimization in AUC with type = "integration"
  ),
  plot = TRUE,
  verbose = FALSE,
  degf = NULL
) {
  stopifnot(
    type %in% c("function", "samples", "grid"),
    utility_type %in% c("AUC", "MSW"),
    is.numeric(nbins),
    is.finite(nbins),
    nbins > 1,
    is.logical(plot),
    is.logical(verbose)
  )
  ## degf deprecation
  df <- .resolve_df(
    df = df,
    degf = degf,
    default = c(1, 5, 20),
    caller = "pseudo_opt"
  )
  if (family == "beta") {
    stop("Optimazation for beta pseudo-targets not implemented.")
  }

  if (is.null(samples)) {
    inits <- c(loc = 0.5, sc = 2.0)
  } else {
    inits <- c(loc = mean(samples), sc = sd(samples))
  }

  get_util <- function(
    pars,
    type,
    log_target,
    samples,
    family,
    df,
    lb,
    ub,
    utility_type,
    nbins,
    verbose
  ) {
    loc <- pars[1]
    sc <- pars[2]

    if (sc <= 0.0) {
      out <- -1.0
    } else {
      pseu <- pseudo_list(
        family = family,
        params = list(loc = loc, sc = sc, df = df),
        lb = lb,
        ub = ub
      )

      if (verbose) {
        cat("trying", pseu$txt, "\n")
      }

      # this function is found in utility_shrinkslice.R; output is the same as utility_shrinkslice()
      out <- utility_pseudo(
        pseudo = pseu,
        log_target = log_target,
        samples = samples,
        type = type,
        nbins = nbins,
        plot = FALSE,
        utility_type = utility_type,
        options = utility_options
      )

      if (verbose) {
        print(out)
        cat("\n")
      }

      if (out > 1.0) {
        # then the result isn't viable. Penalize.
        out <- -1.0
      }
    }

    out
  }

  opt <- list()

  if (family != "t") {
    df <- NA
  }

  for (i in 1:length(df)) {
    opt[[i]] <- optim(
      inits,
      get_util,
      control = list(fnscale = -1, reltol = tol_opt),
      type = type,
      log_target = log_target,
      samples = samples,
      family = family,
      df = df[i],
      lb = lb,
      ub = ub,
      utility_type = utility_type,
      nbins = nbins,
      verbose = verbose
    )
  }

  opt_converged <- sapply(opt, function(obj) obj$convergence == 0)
  if (any(opt_converged)) {
    if (any(!opt_converged)) {
      warning(
        "Pseudo-target utility optimization failed for df = ",
        paste(which(isFALSE(opt_converged)), collapse = " ")
      )
    }
    indx_converged <- which(opt_converged)
    use_indx <- indx_converged[which.max(sapply(indx_converged, function(k) {
      opt[[k]]$value
    }))]
  } else {
    stop("All pseudo-target utility optimizations failed.")
  }

  pseu <- pseudo_list(
    family = family,
    params = list(
      loc = unname(opt[[use_indx]]$par[1]),
      sc = unname(opt[[use_indx]]$par[2]),
      df = df[use_indx]
    ),
    lb = lb,
    ub = ub
  )

  util <- utility_pseudo(
    pseudo = pseu,
    log_target = log_target,
    samples = samples,
    type = type,
    nbins = nbins,
    plot = plot,
    utility_type = utility_type,
    options = utility_options
  )

  list(
    pseudo = pseu,
    utility = util,
    utility_type = utility_type,
    opt = opt[[use_indx]],
    nbins = nbins,
    utility_options = utility_options,
    tol_opt = tol_opt
  )
}
