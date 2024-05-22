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
#'  \code{t}: location \code{loc}, scale \code{sc}, and degrees of freedom \code{degf}
#'
#'  \code{cauchy}: location \code{loc} and scale \code{sc}
#'
#'  \code{norm}: location \code{loc} and scale \code{sc}
#'
#'  \code{logistic}: location \code{loc} and scale \code{sc}
#'
#'  \code{beta}: scale \code{scale1} and scale \code{scale2}
#'
#' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}. Not operative in family \code{beta}.
#' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}. Not operative in family \code{beta}.
#' @param log_p (Not implemented) Logical: evaluate distribution and quantile functions using the log probability.
#' @param name String appending optional message to the textual name of the distribution.
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
#' @references
#' Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," *arXiv preprint arXiv:###*
#'
#' @importFrom stats pt qt dt
#' @importFrom stats pcauchy qcauchy dcauchy
#' @importFrom stats pnorm qnorm dnorm
#' @importFrom stats plogis qlogis dlogis
#' @importFrom stats pbeta qbeta dbeta
#'
#' @export
#' @examples
#' pseu <- pseudo_list(family = "t", params = list(loc = 0.0, sc = 1.0, degf = 5),
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
pseudo_list <- function(family, params, lb = -Inf, ub = Inf,
                        log_p = FALSE, name = NULL) {

  if (family == "t") {

    if (params$degf == 1) {

      out <- pseudo_cauchy_list(loc = params$loc, sc = params$sc,
                                lb = lb, ub = ub, log_p = log_p, name = name)
      out$params <- params

    } else {
      out <- pseudo_t_list(loc = params$loc, sc = params$sc, degf = params$degf,
                           lb = lb, ub = ub, log_p = log_p, name = name)
    }

  } else if (family == "cauchy") {

    out <- pseudo_cauchy_list(loc = params$loc, sc = params$sc,
                              lb = lb, ub = ub, log_p = log_p, name = name)

  } else if (family == "normal") {

    out <- pseudo_normal_list(loc = params$loc, sc = params$sc,
                              lb = lb, ub = ub, log_p = log_p, name = name)

  } else if (family == "logistic") {

    out <- pseudo_logistic_list(loc = params$loc, sc = params$sc,
                                lb = lb, ub = ub, log_p = log_p, name = name)

  } else if (family == "beta") {

    out <- pseudo_beta_list(shape1 = params$shape1, shape2 = params$shape2,
                            log_p = log_p, name = name)

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
# #' @param degf Positive numeric scalar giving the degrees of freedom parameter.
# #' @param lb Numeric scalar giving the value of left truncation. Defaults to \code{-Inf}.
# #' @param ub Numeric scalar giving the value of right truncation. Defaults to \code{Inf}.
# #' @param log_p (Not implemented) Logical: evaluate distribution and quantile functions using the log probability.
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
pseudo_t_list <- function(loc, sc, degf, lb = -Inf, ub = Inf, log_p = FALSE, name = NULL) {

  txt <- paste0("t(loc = ", round(loc,2), ", sc = ", round(sc,2), ", degf = ", round(degf), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  plb <- pt((lb - loc)/sc, df = degf)
  pub <- pt((ub - loc)/sc, df = degf)
  normc <- pub - plb
  lognormc <- log(normc)

  logsc <- log(sc)

  list(
    d = function(x) {
      if (x > lb && x < ub) {
        out <- dt((x - loc) / sc, df = degf) / sc / normc
      } else {
        out <- 0.0
      }
      out
    },
    ld = function(x) {
      if (x > lb && x < ub) {
        out <- dt((x - loc) / sc, df = degf, log = TRUE) - logsc - lognormc
      } else {
        out <- -Inf
      }
      out
    },
    q = function(u, log.p = FALSE) {
      qt(plb + u * normc, log.p = log.p, df = degf) * sc + loc
    },
    p = function(x) {
      if (x < lb) {
        out <- 0.0
      } else if (x <= ub) {
        out <- (pt((x - loc) / sc, df = degf) - plb) / normc
      } else {
        out <- 1.0
      }
      out
    },
    txt = txt,
    params = list(loc = loc, sc = sc, degf = degf),
    lb = lb, ub = ub
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
# #' @param log_p (Not implemented) Logical: evaluate distribution and quantile functions using the log probability.
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
pseudo_cauchy_list <- function(loc, sc, lb = -Inf, ub = Inf, log_p = FALSE, name = NULL) {

  txt <- paste0("Cauchy(loc = ", round(loc,2), ", sc = ", round(sc,2), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  plb <- pcauchy((lb - loc)/sc)
  pub <- pcauchy((ub - loc)/sc)
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
    lb = lb, ub = ub
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
# #' @param log_p (Not implemented) Logical: evaluate distribution and quantile functions using the log probability.
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
pseudo_normal_list <- function(loc, sc, lb = -Inf, ub = Inf, log_p = FALSE, name = NULL) {

  txt <- paste0("normal(loc = ", round(loc,2), ", sc = ", round(sc,2), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  plb <- pnorm((lb - loc)/sc)
  pub <- pnorm((ub - loc)/sc)
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
    params = list( loc = loc, sc = sc),
    lb = lb, ub = ub
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
# #' @param log_p (Not implemented) Logical: evaluate distribution and quantile functions using the log probability.
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
pseudo_logistic_list <- function(loc, sc, lb = -Inf, ub = Inf, log_p = FALSE, name = NULL) {

  txt <- paste0("logistic(loc = ", round(loc,2), ", sc = ", round(sc,2), ")")
  if (!is.null(name)) {
    txt <- paste0(txt, ", ", name)
  }

  if (lb > -Inf || ub < Inf) {
    txt <- paste0(txt, " I(", lb, " < x < ", ub, ")")
  }

  plb <- plogis((lb - loc)/sc)
  pub <- plogis((ub - loc)/sc)
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
    lb = lb, ub = ub
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
# #' @param log_p (Not implemented) Logical: evaluate distribution and quantile functions using the log probability.
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
pseudo_beta_list <- function(shape1, shape2, log_p = FALSE, name = NULL) {

  txt <- paste0("beta(shape1 = ", round(shape1,2), ", shape2 = ", round(shape2,2), ")")
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
    lb = 0.0, ub = 1.0
  )
}
