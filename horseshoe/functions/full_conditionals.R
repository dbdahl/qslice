
upd_lamj2 <- function(lamj2_old, betaj, sig2, tau2) {
  ## betaj / (sig * tau)  ~ N(0, lamj2)
  ## sqrt(lam2) ~ Cauchy+(scale = 1)
  ## update for lamj2
  ## augmentation with eta ~ IG(1/2, 1); lamj2 ~ IG(1/2, 1/eta)

  ssx <- betaj^2 / (sig2 * tau2)

  eta <- 1.0 / rgamma(1.0, shape = 1.0, rate = (1.0 / lamj2_old + 1.0))
  a1 <- 1.0 # since n = 1 in this context
  b1 <- 0.5 * ssx + 1.0 / eta

  1.0 / rgamma(1, shape = a1, rate = b1)
}

constr_Mtau <- function(tau2, lam2, X, n) {

  XD <- X * rep(lam2, each = n)
  XDXt <- tcrossprod(XD, X)

  Mtau <- tau2 * XDXt + diag(n)

  Mtau
}

sig2_fc <- function(y, X, n0, s02, tau2, lam2) {

  n <- length(y)
  a1 <- 0.5 * (n0 + n)

  Mtau <- constr_Mtau(tau2 = tau2, lam2 = lam2, X = X, n = n)

  L <- chol(Mtau) |> t()

  Linvy <- forwardsolve(L, y)
  ss <- crossprod(Linvy) |> drop()
  b1 <- 0.5 * (n0*s02 + ss)

  list(shape = a1, scale = b1, Mtau = Mtau, Linvy = Linvy, ss = ss)
}

upd_sig2 <- function(y, X, n0, s02, tau2, lam2) {
  fc <- sig2_fc(y = y, X = X, n0 = n0, s02 = s02, tau2 = tau2, lam2 = lam2)
  1.0 / rgamma(1, shape = fc$shape, rate = fc$scale)
}

upd_beta <- function(y, X, tau2, lam2, sig2) {

  ## Fast update for regression coefficients in global-local shrinkage by Bhattacharya et al. (2016)

  n <- length(y)
  p <- ncol(X)
  sig <- sqrt(sig2)
  Dvec <- tau2 * lam2

  uu <- rnorm(p) * sqrt(Dvec)
  ff <- rnorm(n)
  vv <- drop(X %*% uu) + ff

  Mtau <- constr_Mtau(tau2 = tau2, lam2 = lam2, X = X, n = n)

  ystar <- y / sig - vv
  vstar <- solve(Mtau, ystar)


  beta0 <- uu + drop(Dvec * crossprod(X, vstar))

  sig * beta0
}


lmarg_tau2 <- function(y, X, tau2, lam2, n0, s02, upper = 1.0e9) {

  n <- length(y)

  if (tau2 < upper) {

    abM <- sig2_fc(y = y, X = X, n0 = n0, s02 = s02, tau2 = tau2, lam2 = lam2)

    # ldetM <- determinant(abM$Mtau$Mtau, logarithm = TRUE)$modulus # can calculate more efficiently?

    ## this alternative calculation actually doesn't help in the n = 100, p = 4000 case... But it does in the n = 500, p = 60 case
    lam <- sqrt(lam2)
    XDh <- X * rep(lam, each = n)
    svdXDh <- svd(XDh)
    indx_nz <- which(svdXDh$d > 0.0)

    eigens_nz <- svdXDh$d[indx_nz]^2
    ldetM <- sum(log1p(tau2 * eigens_nz))

    lfc_tau <- -0.5 * ldetM - abM$shape * log(abM$scale) - log1p(tau2)
    ljacob <- -0.5 * log(tau2) # this is for the transformation tau -> tau2

    out <- lfc_tau + ljacob

  } else {
    warning(paste0("truncated tau2, which was proposed to be ", tau2, ", at ", upper)) # this is necessary for stable computation of target for tau2
    out <- -Inf
  }

  out
}


update_tau2 <- function(state, prior, data, sampler_tau2, upper = 1.0e9) {

  ## state is a list with: beta, sig2, lam2, tau2, and latent_s (latent slice for tau2), and iter
  ## prior is a list with: n0, s20
  ## data is a list with: X, y, n = length(y) = nrow(X), p = ncol(X)

  if (isTRUE(sampler_tau2$logscale)) {

    ltarget <- function(ltau2) {

      lp_tau2 <- lmarg_tau2(y = data$y, X = data$X,
                            tau2 = exp(ltau2),
                            lam2= state$lam2,
                            n0 = prior$n0, s02 = prior$s02,
                            upper = upper)

      lp_tau2 + ltau2 # add Jacobian for transformation tau2 -> log(tau2)
    }

  } else {

    ltarget <- function(tau2) {

      if (tau2 > 0.0) {

        out <- lmarg_tau2(y = data$y, X = data$X,
                          tau2 = tau2, lam2= state$lam2,
                          n0 = prior$n0, s02 = prior$s02,
                          upper = upper)

      } else {
        out <- -Inf
      }

      out
    }

  }

  old_val <- ifelse(sampler_tau2$logscale, log(state$tau2), state$tau2)
  if(isTRUE(sampler_tau2$logscale)) {
    support <- c(-Inf, Inf)
  } else {
    support <- c(0.0, Inf)
  }

  if (sampler_tau2$type == "rw") {

    tmp <- random_walk_sampler(lf = ltarget, support = support,
                               x_0 = old_val, c = sampler_tau2$c)
    tmp$nEvaluations <- 2

  } else if (sampler_tau2$type == "stepping") {

    tmp <- slice_stepping_out(x = old_val,
                              log_target = ltarget,
                              w = sampler_tau2$w)

  } else if (sampler_tau2$type == "gess") {

    tmp <- slice_genelliptical(x = old_val,
                               log_target = ltarget,
                               mu = sampler_tau2$loc,
                               sigma = sampler_tau2$sc,
                               df = sampler_tau2$degf)

  } else if (sampler_tau2$type == "latent") {

    tmp <- slice_latent(x = old_val, s = state$latent_s,
                        log_target = ltarget,
                        rate = sampler_tau2$rate)
    state$latent_s <- tmp$s

  } else if (sampler_tau2$type == "imh") {

    tmp <- imh_pseudo(x = old_val,
                      log_target = ltarget,
                      pseudo = sampler_tau2$pseudo)

    tmp$nEvaluations <- 2

  } else if (sampler_tau2$type == "Qslice") {

    tmp <- slice_quantile(x = old_val,
                          log_target = ltarget,
                          pseudo = sampler_tau2$pseudo)
  }

  state$tau2 <- ifelse(sampler_tau2$logscale, exp(tmp$x), tmp$x)

  list(tau2 = state$tau2, state = state, extras = tmp)
}
