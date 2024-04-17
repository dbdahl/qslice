## Authors: Matt Heiner, Sam Johnson

# source("./MH_samplers.R")

update_beta <- function(state, prior, data) {

  ## state is a list with: g, psi
  ## prior is a list with: beta_0
  ## data is a list with: inv_chol_XtX inv(chol(XtX)),
  ##    beta_mle, p = ncol(X)

  q <- state$g / (1.0 + state$g)
  beta_mean <- q * data$beta_mle + (1.0 - q) * prior$beta_0
  z <- rnorm(data$p)
  out <- sqrt(q / state$psi) * data$inv_chol_XtX %*% z + beta_mean
  drop(out)
}

update_psi <- function(state, prior, data) {

  ## state is a list with: g, beta
  ## prior is a list with: a_0, b_0, beta_0
  ## data is a list with: n, y, X, XtX

  a_n <- prior$a_0 + 0.5*(data$n + data$p)
  residuals <- data$y - drop(data$X %*% state$beta)
  beta_diff <- state$beta - prior$beta_0
  b_n <- prior$b_0 + 0.5 * sum(residuals^2) +
    0.5 * drop(beta_diff %*% data$XtX %*% beta_diff) / state$g
  rgamma(1, a_n, b_n)
}

update_g <- function(state, prior, data, sampler) {

  ## state is list with: g (current value), psi, beta, and latent_s (latent slice)
  ## prior is list with: beta_0, g_max
  ## data is list with: X, logdet_XtX, p = ncol(X)

  qq0 <- data$X %*% (state$beta - prior$beta_0) |> drop()
  qq1 <- drop(crossprod(qq0))

  if (sampler$logG) {
    ltarget <- function(lgg) {
      ## Likelihood for g is unnormalized Gaussian density for beta where
      ##  mean is beta_0 and
      ##  Cov is inv_XtX * g / psi
      ## Prior on g is hyper-g (Liang et al, 2008) on (0, g_max)

      if (lgg <= log(prior$g_max)) {
        qq <- state$psi * qq1 / exp(lgg)
        # logdet <- (data$p * log(gg)) - data$logdet_XtX - (data$p * log(state$psi))
        logdet <- data$p * lgg # only this part is a fn of g
        lpri <- prior$a_g * log1p(exp(lgg))
        out <- -0.5 * (logdet + lpri + qq) + lgg
      } else {
        out <- -Inf
      }

      out
    }
  } else {
    ltarget <- function(gg) {
      ## Likelihood for g is unnormalized Gaussian density for beta where
      ##  mean is beta_0 and
      ##  Cov is inv_XtX * g / psi
      ## Prior on g is hyper-g (Liang et al, 2008) on (0, g_max)

      if ((gg > 0.0) & (gg <= prior$g_max)) {
        qq <- state$psi * qq1 / gg
        # logdet <- (data$p * log(gg)) - data$logdet_XtX - (data$p * log(state$psi))
        logdet <- data$p * log(gg) # only this part is a fn of g
        lpri <- prior$a_g * log(gg + 1.0)
        out <- -0.5 * (logdet + lpri + qq)
      } else {
        out <- -Inf
      }

      out
    }
  }

  g_old <- ifelse(sampler$logG, log(state$g), state$g)
  support <- c(0.0, prior$g_max)
  if(sampler$logG) {
    support <- log(support)
  }

  if (sampler$subtype %in% c("Laplace", "Laplace_wide")) {

    tmp_pseu <- lapproxt(lf = ltarget, init = g_old,
                         sc_adj = sampler$sc_adj,
                         lb = support[1], ub = support[2],
                         maxit = sampler$maxit)

    sampler[["loc"]] <- tmp_pseu$loc
    sampler[["sc"]] <- tmp_pseu$sc
    sampler[["degf"]] <- tmp_pseu$degf

    sampler[["pseudo_lpdf"]] <- tmp_pseu$ld
    sampler[["pseudo_inv_cdf"]] <- tmp_pseu$q

  } else if (sampler$subtype %in% c("Laplace_analytic", "Laplace_analytic_wide")) {

    tmp_pseu <- lapxt_g(p = data$p, psi = state$psi, Q = qq1, a = prior$a_g,
                        sc_adj = sampler$sc_adj, degf = sampler$degf,
                        lb = support[1], ub = support[2],
                        logG = sampler$logG)

    sampler[["loc"]] <- tmp_pseu$loc
    sampler[["sc"]] <- tmp_pseu$sc

    sampler[["pseudo_lpdf"]] <- tmp_pseu$ld
    sampler[["pseudo_inv_cdf"]] <- tmp_pseu$q

  } else if (sampler$subtype %in% c("MM", "MM_wide")) {

    stopifnot(isFALSE(sampler$logG))

    A <- 0.5 * (prior$a_g + data$p - 2)
    B <- 0.5 * state$psi * qq1

    sampler[["loc"]] <- B / (A + 1.0) # IG mode
    sampler[["sc"]] <- sampler$sc_adj * B / ((A - 1.0) * sqrt(A - 2.0))
    sampler[["degf"]] <- sampler$degf

    tmp_pseu <- pseudo_t_list(loc = sampler$loc, sc = sampler$sc, degf = sampler$degf,
                              lb = support[1], ub = support[2])

    sampler[["pseudo_lpdf"]] <- tmp_pseu$ld
    sampler[["pseudo_inv_cdf"]] <- tmp_pseu$q

  }

  if (sampler$type == "Gibbs") {

    stopifnot(isFALSE(sampler$logG) && prior$a_g == 0)

    a1 <- data$p / 2.0 - 1.0
    b1 <- state$psi * qq1 / 2.0

    u_max <- pinvgamma(prior$g_max, shape = a1, scale = b1)
    u <- runif(1, min = 0.0, max = u_max)
    tmp <- list( x = qinvgamma(u, shape = a1, scale = b1),
                 nEvaluations = 0)

  } else if (sampler$type == "rw") {

    tmp <- random_walk_sampler(lf = ltarget, support = support,
                               x_0 = g_old, sampler$c)
    tmp$nEvaluations <- 2

  } else if (sampler$type == "stepping") {

    tmp <- slice_sampler_stepping_out(x = g_old,
                                      target = ltarget,
                                      w = sampler$w)

  } else if (sampler$type == "gess") {

    tmp <- slice_sampler_generalized_elliptical(x = g_old,
                                                target = ltarget,
                                                mu = sampler$loc,
                                                sigma = sampler$sc,
                                                df = sampler$degf)

  } else if (sampler$type == "latent") {

    tmp <- slice_sampler_latent(x = g_old, s = state$latent_s,
                                target = ltarget,
                                rate = sampler$rate)
    state$latent_s <- tmp$s

  } else if (sampler$type == "imh") {

    tmp <- IMH_sampler(lf = ltarget, x_0 = g_old,
                       pseudo_lpdf = sampler$pseudo_lpdf,
                       pseudo_inv_cdf = sampler$pseudo_inv_cdf)

    tmp$nEvaluations <- 2

  } else if (sampler$type == "Qslice") {

    tmp <- slice_sampler_transform(x = g_old,
                                   target = ltarget,
                                   pseudo_log_pdf = sampler$pseudo_lpdf,
                                   pseudo_inv_cdf = sampler$pseudo_inv_cdf)
  }

  state$g <- ifelse(sampler$logG, exp(tmp$x), tmp$x)

  list(g = state$g, state = state, extras = tmp)
}
