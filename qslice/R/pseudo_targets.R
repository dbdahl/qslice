## degf deprecation
.normalize_df_params <- function(params, caller = "pseudo_list") {
  has_df <- !is.null(params$df)
  has_degf <- !is.null(params$degf)

  if (has_df && has_degf) {
    stop(
      caller,
      "(): `params` must contain only `df`; `degf` is deprecated.",
      call. = FALSE
    )
  }

  if (has_degf) {
    warning(
      caller,
      "(): `params$degf` is deprecated; use `params$df`.",
      call. = FALSE
    )

    params$df <- params$degf
    params$degf <- NULL
  }

  params
}

.log1mexp <- function(x) {
  ## Calculates log(1 - exp(x)) for x <= 0.
  stopifnot(all(x <= 0.0))

  out <- numeric(length(x))
  use_log1p <- x < -0.6931472 # -0.6931472 = -log(2)

  out[use_log1p] <- log1p(-exp(x[use_log1p]))
  out[!use_log1p] <- log(-expm1(x[!use_log1p]))

  out
}

.log_addexp <- function(lx, ly) {
  ## Calculates log(exp(lx) + exp(ly)).
  ## From R stats src/nmath/pgamma.c
  m <- pmax(lx, ly)
  out <- m + log1p(exp(-abs(lx - ly)))
  out[is.infinite(m) & m < 0.0] <- -Inf
  out
}

.log_subexp <- function(lx, ly) {
  ## Calculates log(exp(lx) - exp(ly)).
  ## From R stats src/nmath/pgamma.c
  if (any(ly > lx)) {
    stop("`ly` must not exceed `lx`.")
  }
  lx + .log1mexp(ly - lx)
}


.pseudo_locscale_log_probs <- function(lb_std, ub_std, fn_cdf_std) {
  ## Log probabilities to the left of, within, and to the right of
  ## the truncation interval on the standardized scale.
  stopifnot(
    length(lb_std) == 1L,
    length(ub_std) == 1L,
    is.function(fn_cdf_std),
    !is.na(lb_std),
    !is.na(ub_std),
    lb_std < ub_std
  )

  if (lb_std == -Inf) {
    lp_left <- -Inf
  } else if (lb_std <= 0.0) {
    lp_left <- fn_cdf_std(lb_std, lower.tail = TRUE, log.p = TRUE)
  } else {
    lp_left <- .log1mexp(
      fn_cdf_std(lb_std, lower.tail = FALSE, log.p = TRUE)
    )
  }

  if (ub_std == Inf) {
    lp_right <- -Inf
  } else if (ub_std >= 0.0) {
    lp_right <- fn_cdf_std(ub_std, lower.tail = FALSE, log.p = TRUE)
  } else {
    lp_right <- .log1mexp(
      fn_cdf_std(ub_std, lower.tail = TRUE, log.p = TRUE)
    )
  }

  if (ub_std <= 0.0) {
    lp_ub <- fn_cdf_std(ub_std, lower.tail = TRUE, log.p = TRUE)
    lp_interval <- .log_subexp(lp_ub, lp_left)
  } else if (lb_std >= 0.0) {
    lp_lb <- fn_cdf_std(lb_std, lower.tail = FALSE, log.p = TRUE)
    lp_interval <- .log_subexp(lp_lb, lp_right)
  } else {
    lp_interval <- .log1mexp(.log_addexp(lp_left, lp_right))
  }

  list(
    lp_left = lp_left,
    lp_interval = lp_interval,
    lp_right = lp_right
  )
}


.pseudo_locscale_components <- function(
  loc,
  sc,
  lb,
  ub,
  fn_ldens_std,
  fn_cdf_std,
  fn_invcdf_std
) {
  ## Construct a truncated location-scale density, CDF, and quantile from
  ## standardized distribution functions.
  stopifnot(
    length(loc) == 1L,
    length(sc) == 1L,
    length(lb) == 1L,
    length(ub) == 1L,
    is.finite(loc),
    is.finite(sc),
    sc > 0.0,
    !is.na(lb),
    !is.na(ub),
    lb < ub,
    is.function(fn_ldens_std),
    is.function(fn_cdf_std),
    is.function(fn_invcdf_std)
  )

  lb_std <- (lb - loc) / sc
  ub_std <- (ub - loc) / sc
  log_probs <- .pseudo_locscale_log_probs(
    lb_std = lb_std,
    ub_std = ub_std,
    fn_cdf_std = fn_cdf_std
  )

  lp_left <- log_probs$lp_left
  lp_interval <- log_probs$lp_interval
  lp_right <- log_probs$lp_right
  logsc <- log(sc)

  fn_ldens_one <- function(x) {
    if (x > lb && x < ub) {
      fn_ldens_std((x - loc) / sc) - logsc - lp_interval
    } else {
      -Inf
    }
  }

  fn_ldens <- function(x) {
    vapply(x, fn_ldens_one, numeric(1))
  }

  fn_dens <- function(x) exp(fn_ldens(x))

  fn_cdf_one <- function(x) {
    if (x <= lb) {
      return(0.0)
    }
    if (x >= ub) {
      return(1.0)
    }

    x_std <- (x - loc) / sc
    if (x_std <= 0.0) {
      lp_x <- fn_cdf_std(x_std, lower.tail = TRUE, log.p = TRUE)
      lp_numerator <- .log_subexp(lp_x, lp_left)
      out <- exp(lp_numerator - lp_interval)
    } else {
      lp_x <- fn_cdf_std(x_std, lower.tail = FALSE, log.p = TRUE)
      lp_upper <- .log_subexp(lp_x, lp_right) - lp_interval
      out <- -expm1(lp_upper)
    }

    min(1.0, max(0.0, out))
  }

  fn_cdf <- function(x) {
    vapply(x, fn_cdf_one, numeric(1))
  }

  fn_invcdf_one <- function(u) {
    if (is.na(u) || u < 0.0 || u > 1.0) {
      stop("Probabilities in `u` must be in [0, 1].", call. = FALSE)
    }
    if (u == 0.0) {
      return(lb)
    }
    if (u == 1.0) {
      return(ub)
    }

    lp_free <- .log_addexp(lp_left, log(u) + lp_interval)
    if (lp_free <= -0.6931471805599453) {
      # log(0.5)
      x_std <- fn_invcdf_std(
        lp_free,
        lower.tail = TRUE,
        log.p = TRUE
      )
    } else {
      lsp_free <- .log_addexp(
        lp_right,
        log1p(-u) + lp_interval
      )
      x_std <- fn_invcdf_std(
        lsp_free,
        lower.tail = FALSE,
        log.p = TRUE
      )
    }

    x_std * sc + loc
  }

  fn_invcdf <- function(u) {
    if (!is.numeric(u)) {
      stop("`u` must be numeric.", call. = FALSE)
    }
    vapply(u, fn_invcdf_one, numeric(1))
  }

  list(
    fn_dens = fn_dens,
    fn_ldens = fn_ldens,
    fn_cdf = fn_cdf,
    fn_invcdf = fn_invcdf,
    lp_left = lp_left,
    lp_interval = lp_interval,
    lp_right = lp_right
  )
}


.pseudo_bounded_log_probs <- function(lb, ub, fn_cdf_free) {
  ## Log probabilities to the left of, within, and to the right of a
  ## truncation interval. Unlike .pseudo_locscale_log_probs(), the stable
  ## representation is selected by tail probability rather than location.
  stopifnot(
    length(lb) == 1L,
    length(ub) == 1L,
    is.finite(lb),
    is.finite(ub),
    lb < ub,
    is.function(fn_cdf_free)
  )

  log_half <- -0.6931471805599453

  lp_lb <- fn_cdf_free(lb, lower.tail = TRUE, log.p = TRUE)
  lsp_lb <- fn_cdf_free(lb, lower.tail = FALSE, log.p = TRUE)
  lp_ub <- fn_cdf_free(ub, lower.tail = TRUE, log.p = TRUE)
  lsp_ub <- fn_cdf_free(ub, lower.tail = FALSE, log.p = TRUE)

  ## Store log(F(lb)) and log(1 - F(ub)) using the smaller tail directly.
  lp_left <- if (lp_lb <= log_half) lp_lb else .log1mexp(lsp_lb)
  lp_right <- if (lsp_ub <= log_half) lsp_ub else .log1mexp(lp_ub)

  if (lp_ub <= log_half) {
    ## The interval is wholly in the lower-probability half.
    lp_interval <- .log_subexp(lp_ub, lp_left)
  } else if (lsp_lb <= log_half) {
    ## The interval is wholly in the upper-probability half.
    lp_interval <- .log_subexp(lsp_lb, lp_right)
  } else {
    ## The interval contains the median: subtract the two excluded tails.
    lp_interval <- .log1mexp(.log_addexp(lp_left, lp_right))
  }

  if (!is.finite(lp_interval)) {
    stop(
      "The truncation interval has zero probability at machine precision.",
      call. = FALSE
    )
  }

  list(
    lp_left = lp_left,
    lp_interval = lp_interval,
    lp_right = lp_right
  )
}


.pseudo_bounded_components <- function(
  lb,
  ub,
  fn_ldens_free,
  fn_cdf_free,
  fn_invcdf_free
) {
  ## Construct a truncated bounded density, CDF, and quantile from the
  ## corresponding untruncated distribution functions.
  stopifnot(
    length(lb) == 1L,
    length(ub) == 1L,
    is.finite(lb),
    is.finite(ub),
    lb < ub,
    is.function(fn_ldens_free),
    is.function(fn_cdf_free),
    is.function(fn_invcdf_free)
  )

  log_probs <- .pseudo_bounded_log_probs(
    lb = lb,
    ub = ub,
    fn_cdf_free = fn_cdf_free
  )

  lp_left <- log_probs$lp_left
  lp_interval <- log_probs$lp_interval
  lp_right <- log_probs$lp_right
  log_half <- -0.6931471805599453

  fn_ldens_one <- function(x) {
    if (is.na(x)) {
      return(x)
    }
    if (x > lb && x < ub) {
      fn_ldens_free(x) - lp_interval
    } else {
      -Inf
    }
  }

  fn_ldens <- function(x) {
    vapply(x, fn_ldens_one, numeric(1))
  }

  fn_dens <- function(x) exp(fn_ldens(x))

  fn_cdf_one <- function(x) {
    if (is.na(x)) {
      return(x)
    }
    if (x <= lb) {
      return(0.0)
    }
    if (x >= ub) {
      return(1.0)
    }

    lp_x <- fn_cdf_free(x, lower.tail = TRUE, log.p = TRUE)
    if (lp_x <= log_half) {
      lp_numerator <- .log_subexp(lp_x, lp_left)
      out <- exp(lp_numerator - lp_interval)
    } else {
      lsp_x <- fn_cdf_free(x, lower.tail = FALSE, log.p = TRUE)
      lp_upper <- .log_subexp(lsp_x, lp_right) - lp_interval
      out <- -expm1(lp_upper)
    }

    min(1.0, max(0.0, out))
  }

  fn_cdf <- function(x) {
    vapply(x, fn_cdf_one, numeric(1))
  }

  fn_invcdf_one <- function(u, log.p) {
    if (is.na(u)) {
      stop("`u` must not be NA or NaN.", call. = FALSE)
    }

    if (log.p) {
      if (u > 0.0) {
        stop("Log probabilities in `u` must not exceed zero.", call. = FALSE)
      }
      if (u == -Inf) {
        return(lb)
      }
      if (u == 0.0) {
        return(ub)
      }
      lp_u <- u
      lsp_u <- .log1mexp(u)
    } else {
      if (u < 0.0 || u > 1.0) {
        stop("Probabilities in `u` must be in [0, 1].", call. = FALSE)
      }
      if (u == 0.0) {
        return(lb)
      }
      if (u == 1.0) {
        return(ub)
      }
      lp_u <- log(u)
      lsp_u <- log1p(-u)
    }

    lp_free <- .log_addexp(lp_left, lp_u + lp_interval)
    if (lp_free <= log_half) {
      fn_invcdf_free(lp_free, lower.tail = TRUE, log.p = TRUE)
    } else {
      lsp_free <- .log_addexp(lp_right, lsp_u + lp_interval)
      fn_invcdf_free(lsp_free, lower.tail = FALSE, log.p = TRUE)
    }
  }

  fn_invcdf <- function(u, log.p = FALSE) {
    if (
      !is.numeric(u) ||
        !is.logical(log.p) ||
        length(log.p) != 1L ||
        is.na(log.p)
    ) {
      stop(
        "`u` must be numeric and `log.p` must be TRUE or FALSE.",
        call. = FALSE
      )
    }
    vapply(u, fn_invcdf_one, numeric(1), log.p = log.p)
  }

  list(
    fn_dens = fn_dens,
    fn_ldens = fn_ldens,
    fn_cdf = fn_cdf,
    fn_invcdf = fn_invcdf,
    lp_left = lp_left,
    lp_interval = lp_interval,
    lp_right = lp_right
  )
}


#' Specify a pseudo-target within a given class
#'
#' Create a list of functions to evaluate a pseudo-target in a given class
#' with supplied parameters (usually location and scale). The distribution is optionally
#' truncated to specified bounds (and renormalized). See Heiner et al. (2026+).
#'
#' The supported classes of pseudo-targets include: \code{t}, \code{cauchy},
#' \code{normal}, \code{logistic}, and \code{beta}.
#'
#' @param family String identifying the distribution family. One of \code{t}, \code{cauchy},
#' \code{normal}, \code{logistic}, and \code{beta}.
#' @param params Named list identifying parameters, which vary by distribution family.
#'
#'  \code{t}: location \code{loc}, scale \code{sc}, and degrees of freedom \code{df}
#'
#'  \code{cauchy}: location \code{loc} and scale \code{sc}
#'
#'  \code{norm}: location \code{loc} and scale \code{sc}
#'
#'  \code{logistic}: location \code{loc} and scale \code{sc}
#'
#'  \code{beta}: shape \code{shape1} and shape \code{shape2}
#'
#' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
#' For family \code{beta}, this is intersected with its natural lower bound of zero.
#' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
#' For family \code{beta}, this is intersected with its natural upper bound of one.
#' @param name String appending optional message to the textual name of the distribution.
#' @returns A list with named components:
#'
#'  \code{d}: function to evaluate the density (finite boundary points of the support evaluate to \code{0})
#'
#'  \code{ld}: function to evaluate the log density (finite boundary points of the support evaluate to \code{-Inf})
#'
#'  \code{q}: function to evaluate the quantile function
#'
#'  \code{p}: function to evaluate the distribution function
#'
#'  \code{txt}: text description of the distribution
#'
#'  \code{params}: repeats the \code{params} argument
#'
#'  \code{lb}: lower boundary of support
#'
#'  \code{ub}: upper boundary of support
#'
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2026+), "Quantile Slice Sampling," *arXiv preprint arXiv:2407.12608* \doi{https://doi.org/10.48550/arXiv.2407.12608}
#'
#' @importFrom stats pt qt dt
#' @importFrom stats pcauchy qcauchy dcauchy
#' @importFrom stats pnorm qnorm dnorm
#' @importFrom stats plogis qlogis dlogis
#' @importFrom stats pbeta qbeta dbeta
#'
#' @export
#' @examples
#' pseu <- pseudo_list(family = "t", params = list(loc = 0.0, sc = 1.0, df = 5),
#'                     lb = 0.0, ub = Inf) # half t
#' str(pseu)
#' pseu$d(1.5)
#' pseu$ld(1.5)
#' pseu$p(1.5)
#' # should match (F(x) - F(0)) / (1 - F(0)) below
#' (pt(1.5, df = 5) - pt(0.0, df = 5)) / pt(0.0, df = 5, lower.tail = FALSE)
#' pseu$q(0.8060963)
#' curve(pseu$d(x), from = -1.0, to = 5.0, n = 1000)
#' pseu <- pseudo_list(family = "cauchy", params = list(loc = 0.0, sc = 1.0),
#'                     lb = 0.0, ub = Inf) # half Cauchy
#' str(pseu)
#' pseu$d(1.5)
#' pseu$ld(1.5)
#' pseu$p(1.5)
#' pseu$q(0.6256659)
#' pseu <- pseudo_list(family = "normal", params = list(loc = 0.0, sc = 1.0),
#'                     lb = 0.0, ub = Inf) # half normal
#' str(pseu)
#' pseu$d(1.5)
#' pseu$ld(1.5)
#' pseu$p(1.5)
#' pseu$q(0.8663856)
#' pseu <- pseudo_list(family = "logistic", params = list(loc = 0.0, sc = 1.0),
#'                     lb = 0.0, ub = Inf) # half logistic
#' str(pseu)
#' pseu$d(1.5)
#' pseu$ld(1.5)
#' pseu$p(1.5)
#' pseu$q(0.635149)
#' pseu <- pseudo_list(family = "beta", params = list(shape1 = 3.0, shape2 = 2.0), lb = 0.2)
#' str(pseu)
#' pseu$d(0.5)
#' pseu$ld(0.5)
#' pseu$p(0.5)
#' # should match (F(x) - F(0.2)) / (1 - F(0.2)) below
#' (pbeta(0.5, 3.0, 2.0) - pbeta(0.2, 3.0, 2.0)) / pbeta(0.2, 3.0, 2.0, lower.tail = FALSE)
#' pseu$q(0.2932771)
#' curve(pseu$d(x), from = 0.0, to = 1.0, n = 1000)
pseudo_list <- function(
  family,
  params,
  lb = -Inf,
  ub = Inf,
  name = NULL
) {
  params <- .normalize_df_params(params) # degf deprecation

  if (family == "t") {
    if (params$df == 1) {
      out <- .pseudo_cauchy_list(
        loc = params$loc,
        sc = params$sc,
        lb = lb,
        ub = ub,
        name = name
      )
      out$params <- params
    } else {
      out <- .pseudo_t_list(
        loc = params$loc,
        sc = params$sc,
        df = params$df,
        lb = lb,
        ub = ub,
        name = name
      )
    }
  } else if (family == "cauchy") {
    out <- .pseudo_cauchy_list(
      loc = params$loc,
      sc = params$sc,
      lb = lb,
      ub = ub,
      name = name
    )
  } else if (family == "normal") {
    out <- .pseudo_normal_list(
      loc = params$loc,
      sc = params$sc,
      lb = lb,
      ub = ub,
      name = name
    )
  } else if (family == "logistic") {
    out <- .pseudo_logistic_list(
      loc = params$loc,
      sc = params$sc,
      lb = lb,
      ub = ub,
      name = name
    )
  } else if (family == "beta") {
    out <- .pseudo_beta_list(
      shape1 = params$shape1,
      shape2 = params$shape2,
      lb = lb,
      ub = ub,
      name = name
    )
  } else {
    stop("Pseudo-target family supplied to pseudo_list() is not supported.")
  }

  out$family <- family
  out
}


# #' Specify a Pseudo-Target within the Student-t Class
# #'
# #' Create a list of functions to evaluate a pseudo-target in the Student-t class
# #' with supplied location, scale, and degrees of freedom. The distribution is optionally
# #' truncated to specified bounds (and renormalized).
# #'
# #'
# #' @param loc Numeric scalar giving the location parameter.
# #' @param sc Positive numeric scalar giving the scale parameter.
# #' @param df Positive numeric scalar giving the degrees of freedom parameter.
# #' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
# #' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
# #' @param name String appending optional message to the textual name of the distribution.
# #' @returns A list with named components:
# #'
# #'  \code{d}: function to evaluate the density (finite boundary points of restricted support evaluate to \code{0})
# #'
# #'  \code{ld}: function to evaluate the log density (finite boundary points of restricted support evaluate to \code{-Inf})
# #'
# #'  \code{q}: function to evaluate the quantile function
# #'
# #'  \code{p}: function to evaluate the distribution function
# #'
# #'  \code{txt}: text description of the distribution
# #'
# #'  \code{params}: returns the parameters passed to the function
# #'
# #'  \code{lb}: lower boundary of support
# #'
# #'  \code{ub}: upper boundary of support
# #'
# #' @importFrom stats pt qt dt
# #' @keywords internal
# #'
.pseudo_t_list <- function(
  loc,
  sc,
  df,
  lb = -Inf,
  ub = Inf,
  name = NULL
) {
  stopifnot(
    length(loc) == 1L,
    length(sc) == 1L,
    length(df) == 1L,
    length(lb) == 1L,
    length(ub) == 1L,
    is.finite(loc),
    is.finite(sc),
    df > 0,
    lb < ub
  )

  txt <- paste0(
    "t(loc = ",
    round(loc, 2),
    ", sc = ",
    round(sc, 2),
    ", df = ",
    round(df),
    ")"
  )
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  fn_ldens_std <- function(x) dt(x, df = df, log = TRUE)
  fn_cdf_std <- function(x, lower.tail = TRUE, log.p = FALSE) {
    pt(x, df = df, lower.tail = lower.tail, log.p = log.p)
  }
  fn_invcdf_std <- function(u, lower.tail = TRUE, log.p = FALSE) {
    qt(u, df = df, lower.tail = lower.tail, log.p = log.p)
  }

  components <- .pseudo_locscale_components(
    loc = loc,
    sc = sc,
    lb = lb,
    ub = ub,
    fn_ldens_std = fn_ldens_std,
    fn_cdf_std = fn_cdf_std,
    fn_invcdf_std = fn_invcdf_std
  )

  list(
    d = components$fn_dens,
    ld = components$fn_ldens,
    q = components$fn_invcdf,
    p = components$fn_cdf,
    txt = txt,
    params = list(loc = loc, sc = sc, df = df),
    lb = lb,
    ub = ub
  )
}


# #' Specify a Cauchy Pseudo-Target
# #'
# #' Create a list of functions to evaluate a Cauchy (Student-t with one degree of freedom) pseudo-target
# #' with supplied location, scale. The distribution is optionally
# #' truncated to specified bounds (and renormalized).
# #'
# #'
# #' @param loc Numeric scalar giving the location parameter.
# #' @param sc Positive numeric scalar giving the scale parameter.
# #' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
# #' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
# #' @param name String appending optional message to the textual name of the distribution.
# #' @returns A list with named components:
# #'
# #'  \code{d}: function to evaluate the density (finite boundary points of restricted support evaluate to \code{0})
# #'
# #'  \code{ld}: function to evaluate the log density (finite boundary points of restricted support evaluate to \code{-Inf})
# #'
# #'  \code{q}: function to evaluate the quantile function
# #'
# #'  \code{p}: function to evaluate the distribution function
# #'
# #'  \code{txt}: text description of the distribution
# #'
# #'  \code{params}: returns the parameters passed to the function
# #'
# #'  \code{lb}: lower boundary of support
# #'
# #'  \code{ub}: upper boundary of support
# #'
# #' @importFrom stats pcauchy qcauchy dcauchy
# #' @keywords internal
# #'
.pseudo_cauchy_list <- function(
  loc,
  sc,
  lb = -Inf,
  ub = Inf,
  name = NULL
) {
  txt <- paste0("Cauchy(loc = ", round(loc, 2), ", sc = ", round(sc, 2), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  fn_ldens_std <- function(x) dcauchy(x, log = TRUE)
  fn_cdf_std <- function(x, lower.tail = TRUE, log.p = FALSE) {
    pcauchy(x, lower.tail = lower.tail, log.p = log.p)
  }
  fn_invcdf_std <- function(u, lower.tail = TRUE, log.p = FALSE) {
    qcauchy(u, lower.tail = lower.tail, log.p = log.p)
  }

  components <- .pseudo_locscale_components(
    loc = loc,
    sc = sc,
    lb = lb,
    ub = ub,
    fn_ldens_std = fn_ldens_std,
    fn_cdf_std = fn_cdf_std,
    fn_invcdf_std = fn_invcdf_std
  )

  list(
    d = components$fn_dens,
    ld = components$fn_ldens,
    q = components$fn_invcdf,
    p = components$fn_cdf,
    txt = txt,
    params = list(loc = loc, sc = sc),
    lb = lb,
    ub = ub
  )
}

# #' Specify a normal Pseudo-Target
# #'
# #' Create a list of functions to evaluate a normal pseudo-target
# #' with supplied location, scale. The distribution is optionally
# #' truncated to specified bounds (and renormalized).
# #'
# #'
# #' @param loc Numeric scalar giving the location parameter.
# #' @param sc Positive numeric scalar giving the scale parameter.
# #' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
# #' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
# #' @param name String appending optional message to the textual name of the distribution.
# #' @returns A list with named components:
# #'
# #'  \code{d}: function to evaluate the density (finite boundary points of restricted support evaluate to \code{0})
# #'
# #'  \code{ld}: function to evaluate the log density (finite boundary points of restricted support evaluate to \code{-Inf})
# #'
# #'  \code{q}: function to evaluate the quantile function
# #'
# #'  \code{p}: function to evaluate the distribution function
# #'
# #'  \code{txt}: text description of the distribution
# #'
# #'  \code{params}: returns the parameters passed to the function
# #'
# #'  \code{lb}: lower boundary of support
# #'
# #'  \code{ub}: upper boundary of support
# #'
# #' @importFrom stats pnorm qnorm dnorm
# #' @keywords internal
# #'
.pseudo_normal_list <- function(
  loc,
  sc,
  lb = -Inf,
  ub = Inf,
  name = NULL
) {
  txt <- paste0("normal(loc = ", round(loc, 2), ", sc = ", round(sc, 2), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  fn_ldens_std <- function(x) dnorm(x, log = TRUE)
  fn_cdf_std <- function(x, lower.tail = TRUE, log.p = FALSE) {
    pnorm(x, lower.tail = lower.tail, log.p = log.p)
  }
  fn_invcdf_std <- function(u, lower.tail = TRUE, log.p = FALSE) {
    qnorm(u, lower.tail = lower.tail, log.p = log.p)
  }

  components <- .pseudo_locscale_components(
    loc = loc,
    sc = sc,
    lb = lb,
    ub = ub,
    fn_ldens_std = fn_ldens_std,
    fn_cdf_std = fn_cdf_std,
    fn_invcdf_std = fn_invcdf_std
  )

  list(
    d = components$fn_dens,
    ld = components$fn_ldens,
    q = components$fn_invcdf,
    p = components$fn_cdf,
    txt = txt,
    params = list(loc = loc, sc = sc),
    lb = lb,
    ub = ub
  )
}


# #' Specify a logistic Pseudo-Target
# #'
# #' Create a list of functions to evaluate a logistic pseudo-target
# #' with supplied location, scale. The distribution is optionally
# #' truncated to specified bounds (and renormalized).
# #'
# #'
# #' @param loc Numeric scalar giving the location parameter.
# #' @param sc Positive numeric scalar giving the scale parameter.
# #' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
# #' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
# #' @param name String appending optional message to the textual name of the distribution.
# #' @returns A list with named components:
# #'
# #'  \code{d}: function to evaluate the density (finite boundary points of restricted support evaluate to \code{0})
# #'
# #'  \code{ld}: function to evaluate the log density (finite boundary points of restricted support evaluate to \code{-Inf})
# #'
# #'  \code{q}: function to evaluate the quantile function
# #'
# #'  \code{p}: function to evaluate the distribution function
# #'
# #'  \code{txt}: text description of the distribution
# #'
# #'  \code{params}: returns the parameters passed to the function
# #'
# #'  \code{lb}: lower boundary of support
# #'
# #'  \code{ub}: upper boundary of support
# #'
# #' @importFrom stats plogis qlogis dlogis
# #' @keywords internal
# #'
.pseudo_logistic_list <- function(
  loc,
  sc,
  lb = -Inf,
  ub = Inf,
  name = NULL
) {
  txt <- paste0("logistic(loc = ", round(loc, 2), ", sc = ", round(sc, 2), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  fn_ldens_std <- function(x) dlogis(x, log = TRUE)
  fn_cdf_std <- function(x, lower.tail = TRUE, log.p = FALSE) {
    plogis(x, lower.tail = lower.tail, log.p = log.p)
  }
  fn_invcdf_std <- function(u, lower.tail = TRUE, log.p = FALSE) {
    qlogis(u, lower.tail = lower.tail, log.p = log.p)
  }

  components <- .pseudo_locscale_components(
    loc = loc,
    sc = sc,
    lb = lb,
    ub = ub,
    fn_ldens_std = fn_ldens_std,
    fn_cdf_std = fn_cdf_std,
    fn_invcdf_std = fn_invcdf_std
  )

  list(
    d = components$fn_dens,
    ld = components$fn_ldens,
    q = components$fn_invcdf,
    p = components$fn_cdf,
    txt = txt,
    params = list(loc = loc, sc = sc),
    lb = lb,
    ub = ub
  )
}


# #' Specify a beta Pseudo-Target
# #'
# #' Create a list of functions to evaluate a beta pseudo-target
# #' with supplied shape parameters.
# #'
# #'
# #' @param shape1 Positive numeric scalar giving the first shape parameter.
# #' @param shape2 Positive numeric scalar giving the second shape parameter.
# #' @param lb Numeric scalar giving the value of left truncation. It is
# #' intersected with the beta distribution's natural lower bound of zero.
# #' @param ub Numeric scalar giving the value of right truncation. It is
# #' intersected with the beta distribution's natural upper bound of one.
# #' @param name String appending optional message to the textual name of the distribution.
# #' @returns A list with named components:
# #'
# #'  \code{d}: function to evaluate the density (finite boundary points of the support evaluate to \code{0})
# #'
# #'  \code{ld}: function to evaluate the log density (finite boundary points of the support evaluate to \code{-Inf})
# #'
# #'  \code{q}: function to evaluate the quantile function
# #'
# #'  \code{p}: function to evaluate the distribution function
# #'
# #'  \code{t}: text description of the distribution
# #'
# #'  \code{params}: returns the parameters passed to the function
# #'
# #' @importFrom stats pbeta qbeta dbeta
# #' @keywords internal
# #'
.pseudo_beta_list <- function(
  shape1,
  shape2,
  lb = -Inf,
  ub = Inf,
  name = NULL
) {
  stopifnot(
    length(shape1) == 1L,
    length(shape2) == 1L,
    length(lb) == 1L,
    length(ub) == 1L,
    is.finite(shape1),
    is.finite(shape2),
    shape1 > 0.0,
    shape2 > 0.0,
    !is.na(lb),
    !is.na(ub),
    lb < ub
  )

  ## Intersect the requested truncation interval with beta's natural support.
  lb <- max(lb, 0.0)
  ub <- min(ub, 1.0)
  if (lb >= ub) {
    stop(
      "The truncation interval must overlap the beta support (0, 1).",
      call. = FALSE
    )
  }

  txt <- paste0(
    "beta(shape1 = ",
    round(shape1, 2),
    ", shape2 = ",
    round(shape2, 2),
    ")"
  )
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > 0.0 || ub < 1.0) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  fn_ldens_free <- function(x) {
    dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)
  }
  fn_cdf_free <- function(x, lower.tail = TRUE, log.p = FALSE) {
    pbeta(
      x,
      shape1 = shape1,
      shape2 = shape2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }
  fn_invcdf_free <- function(u, lower.tail = TRUE, log.p = FALSE) {
    qbeta(
      u,
      shape1 = shape1,
      shape2 = shape2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  components <- .pseudo_bounded_components(
    lb = lb,
    ub = ub,
    fn_ldens_free = fn_ldens_free,
    fn_cdf_free = fn_cdf_free,
    fn_invcdf_free = fn_invcdf_free
  )

  list(
    d = components$fn_dens,
    ld = components$fn_ldens,
    q = components$fn_invcdf,
    p = components$fn_cdf,
    txt = txt,
    params = list(shape1 = shape1, shape2 = shape2),
    lb = lb,
    ub = ub
  )
}
