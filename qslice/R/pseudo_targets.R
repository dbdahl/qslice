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

.log_xplusy <- function(lx, ly) {
  ## Calculates log(exp(lx) + exp(ly))
  ## From R stats src/nmath/pgamma.c
  pmax(lx, ly) + log1p(exp(-abs(lx - ly)))
}

.log_xminusy <- function(lx, ly) {
  ## Calculates log(exp(lx) - exp(ly)).
  ## From R stats src/nmath/pgamma.c
  if (any(ly > lx)) {
    stop("`ly` must not exceed `lx`.")
  }
  lx + .log1mexp(ly - lx)
}


#' Specify a pseudo-target within a given class
#'
#' Create a list of functions to evaluate a pseudo-target in a given class
#' with supplied parameters (usually location and scale). The distribution is optionally
#' truncated to specified bounds (and renormalized). See Heiner et al. (2024+).
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
#' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}. Not operative in family \code{beta}.
#' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}. Not operative in family \code{beta}.
#' @param name String appending optional message to the textual name of the distribution.
#' @param log_tails Logical: evaluate distribution and quantile functions using log probabilities to improve numerical accuracy.
#' The scale of the original variable and its CDF-transform are unchanged. Defaults to true.
#' @returns A list with named components:
#'
#'  \code{d}: function to evaluate the density
#'
#'  \code{ld}: function to evaluate the log density
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
#'  \code{log_tails}: logical variable indicating whether the \code{q} and \code{p} functions support ### I MIGHT BE ABLE TO GET RID OF ALL OF THIS
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
#' pseu$q(0.8060963)
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
#' pseu <- pseudo_list(family = "beta", params = list(shape1 = 2.0, shape2 = 1.0))
#' str(pseu)
#' pseu$d(0.5)
#' pseu$ld(0.5)
#' pseu$p(0.5)
#' pseu$q(0.25)
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
# #'  \code{d}: function to evaluate the density
# #'
# #'  \code{ld}: function to evaluate the log density
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

  if (lb == -Inf) {
    lp_left <- -Inf
    if (ub == Inf) {
      lp_interval <- 0.0
      lp_right <- -Inf
    } else if (ub <= loc) {
      lcdf_ub <- pt((ub - loc) / sc, df = df, lower.tail = TRUE, log.p = TRUE)
      lp_interval <- lcdf_ub
      lp_right <- .log1mexp(lcdf_ub)
    } else if (ub > loc) {
      lsurv_ub <- pt((ub - loc) / sc, df = df, lower.tail = FALSE, log.p = TRUE)
      lp_right <- lsurv_ub
      lp_interval <- .log1mexp(lsurv_ub)
    } else {
      stop("Invalid truncation bounds lb and ub")
    }
  } else if (lb <= loc) {
    lcdf_lb <- pt((lb - loc) / sc, df = df, lower.tail = TRUE, log.p = TRUE)
    lp_left <- lcdf_lb
    if (ub == Inf) {
      lp_interval <- .log1mexp(lcdf_lb)
      lp_right <- -Inf
    } else if (ub <= loc) {
      lcdf_ub <- pt((ub - loc) / sc, df = df, lower.tail = TRUE, log.p = TRUE)
      lp_right <- .log1mexp(lcdf_ub)
      lp_interval <- .log_xminusy(lcdf_ub, lcdf_lb)
    } else if (ub > loc) {
      lsurv_ub <- pt((ub - loc) / sc, df = df, lower.tail = FALSE, log.p = TRUE)
      lp_interval <- .log1mexp(.log_xplusy(lcdf_lb, lsurv_ub))
      lp_right <- lsurv_ub
    } else {
      stop("Invalid truncation bounds lb and ub")
    }
  } else if (lb > loc) {
    lsurv_lb <- pt((lb - loc) / sc, df = df, lower.tail = FALSE, log.p = TRUE)
    lp_left <- .log1mexp(lsurv_lb)
    if (ub == Inf) {
      lp_interval <- lsurv_lb
      lp_right <- -Inf
    } else if (ub <= loc) {
      stop("Invalid truncation bounds lb and ub")
    } else if (ub > loc) {
      lsurv_ub <- pt((ub - loc) / sc, df = df, lower.tail = FALSE, log.p = TRUE)
      lp_interval <- .log_xminusy(lsurv_lb, lsurv_ub)
      lp_right <- lsurv_ub
    } else {
      stop("Invalid truncation bounds lb and ub")
    }
  } else {
    stop("Invalid truncation bounds lb and ub")
  }

  logsc <- log(sc)

  fn_ldens <- function(x) {
    if (x > lb && x < ub) {
      out <- dt((x - loc) / sc, df = df, log = TRUE) - logsc - lp_interval
    } else {
      out <- -Inf
    }
    out
  }

  fn_dens <- function(x) exp(fn_ldens(x))

  fn_cdf <- function(x) {
    if (x < lb) {
      out <- 0.0
    } else if (x <= ub) {
      lcdf_free <- pt(
        (x - loc) / sc,
        df = df,
        lower.tail = TRUE,
        log.p = TRUE
      )
      lnumerator <- .log_xminusy(lcdf_free, lp_left)
      out <- exp(lnumerator - lp_interval)
    } else {
      out <- 1.0
    }
    out
  }

  fn_invcdf <- function(u) {
    logp_free <- .log_xplusy(lp_left, log(u) + lp_interval)

    if (logp_free <= -0.6931472) {
      # -0.6931472 = log(0.5)
      z_out <- qt(logp_free, df = df, lower.tail = TRUE, log.p = TRUE)
    } else {
      log1mp_free <- .log_xplusy(lp_right, log1p(-u) + lp_interval)
      z_out <- qt(log1mp_free, df = df, lower.tail = FALSE, log.p = TRUE)
    }

    z_out * sc + loc
  }

  list(
    d = fn_dens,
    ld = fn_ldens,
    q = fn_invcdf,
    p = fn_cdf,
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
# #'  \code{d}: function to evaluate the density
# #'
# #'  \code{ld}: function to evaluate the log density
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

  plb <- pcauchy((lb - loc) / sc)
  pub <- pcauchy((ub - loc) / sc)
  normc <- pub - plb
  lognormc <- log(normc)

  logsc <- log(sc)

  list(
    d = function(x) {
      if (x > lb && x < ub) {
        out <- dcauchy((x - loc) / sc) / sc / normc
      } else {
        out <- 0.0
      }
      out
    },
    ld = function(x) {
      if (x > lb && x < ub) {
        out <- dcauchy((x - loc) / sc, log = TRUE) - logsc - lognormc
      } else {
        out <- -Inf
      }
      out
    },
    q = function(u, log.p = FALSE) {
      qcauchy(plb + u * normc, log.p = log.p) * sc + loc
    },
    p = function(x) {
      if (x < lb) {
        out <- 0.0
      } else if (x <= ub) {
        out <- (pcauchy((x - loc) / sc) - plb) / normc
      } else {
        out <- 1.0
      }
      out
    },
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
# #'  \code{d}: function to evaluate the density
# #'
# #'  \code{ld}: function to evaluate the log density
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

  plb <- pnorm((lb - loc) / sc)
  pub <- pnorm((ub - loc) / sc)
  normc <- pub - plb
  lognormc <- log(normc)

  logsc <- log(sc)

  list(
    d = function(x) {
      if (x > lb && x < ub) {
        out <- dnorm((x - loc) / sc) / sc / normc
      } else {
        out <- 0.0
      }
      out
    },
    ld = function(x) {
      if (x > lb && x < ub) {
        out <- dnorm((x - loc) / sc, log = TRUE) - logsc - lognormc
      } else {
        out <- -Inf
      }
      out
    },
    q = function(u, log.p = FALSE) {
      qnorm(plb + u * normc, log.p = log.p) * sc + loc
    },
    p = function(x) {
      if (x < lb) {
        out <- 0.0
      } else if (x <= ub) {
        out <- (pnorm((x - loc) / sc) - plb) / normc
      } else {
        out <- 1.0
      }
      out
    },
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
# #'  \code{d}: function to evaluate the density
# #'
# #'  \code{ld}: function to evaluate the log density
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

  plb <- plogis((lb - loc) / sc)
  pub <- plogis((ub - loc) / sc)
  normc <- pub - plb
  lognormc <- log(normc)

  logsc <- log(sc)

  list(
    d = function(x) {
      if (x > lb && x < ub) {
        out <- dlogis((x - loc) / sc) / sc / normc
      } else {
        out <- 0.0
      }
      out
    },
    ld = function(x) {
      if (x > lb && x < ub) {
        out <- dlogis((x - loc) / sc, log = TRUE) - logsc - lognormc
      } else {
        out <- -Inf
      }
      out
    },
    q = function(u, log.p = FALSE) {
      qlogis(plb + u * normc, log.p = log.p) * sc + loc
    },
    p = function(x) {
      if (x < lb) {
        out <- 0.0
      } else if (x <= ub) {
        out <- (plogis((x - loc) / sc) - plb) / normc
      } else {
        out <- 1.0
      }
      out
    },
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
# #' @param name String appending optional message to the textual name of the distribution.
# #' @returns A list with named components:
# #'
# #'  \code{d}: function to evaluate the density
# #'
# #'  \code{ld}: function to evaluate the log density
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
.pseudo_beta_list <- function(shape1, shape2, name = NULL) {
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

  list(
    d = function(x) {
      dbeta(x, shape1 = shape1, shape2 = shape2, log = FALSE)
    },
    ld = function(x) {
      dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)
    },
    q = function(u, log.p = FALSE) {
      qbeta(u, shape1 = shape1, shape2 = shape2, log.p = log.p)
    },
    p = function(x) {
      pbeta(x, shape1 = shape1, shape2 = shape2)
    },
    txt = txt,
    params = list(shape1 = shape1, shape2 = shape2),
    lb = 0.0,
    ub = 1.0
  )
}
