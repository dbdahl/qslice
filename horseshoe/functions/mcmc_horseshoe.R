
mcmc_hs <- function(state, prior, data, sampler, n_iter,
                    n_thin = 1, save = TRUE, prog = 0,
                    upper_tau2 = 1.0e9) {

  ## state is a list with: beta, sig2, lam2, tau2, and latent_s (latent slice for tau2), and iter
  ## prior is a list with: n0, s20
  ## data is a list with: X, y, n = length(y) = nrow(X), p = ncol(X)

  if (save) {
    sims <- rep(list(state), n_iter)
    extras <- rep(list(), n_iter)
  } else {
    sims <- NULL
    extras <- NULL
  }

  for (i in 1:n_iter) {

    for (ii in 1:n_thin) {

      for (j in 1:data$p) {
        state$lam2[j] <- upd_lamj2(lamj2_old = state$lam2[j], betaj = state$beta[j], sig2 = state$sig2, tau2 = state$tau2)
      }

      tmp <- update_tau2(state = state, prior = prior, data = data,
                         sampler = sampler[["tau2"]],
                         upper = upper_tau2)
      state <- tmp$state # includes update of s in latent sampler

      state$sig2 <- upd_sig2(y = data$y, X = data$X, n0 = prior$n0, s02 = prior$s02,
                             tau2 = state$tau2, lam2 = state$lam2)

      state$beta <- upd_beta(y = data$y, X = data$X,
                             tau2 = state$tau2, lam2 = state$lam2, sig2 = state$sig2)


      state$iter <- state$iter + 1

    }

    if (save) {
      sims[[i]] <- state
      extras[[i]] <- tmp[["extras"]]
    }

    if (prog > 0) {
      if (i %% prog == 0) {
        cat("Sample ", i, " of ", n_iter, "\n")
        timestamp()
      }
    }

  }

  list(state = state, prior = prior, data = data, sampler = sampler,
       save = save, sims = sims, extras = extras, n_iter = n_iter)
}


time_hs <- function(state, prior, data, sampler, n_iter, param,
                    upper_tau2 = 1.0e9) {

  require("coda")

  time_out <- system.time({
    mcmc_out <- mcmc_hs(state = state,
                        prior = prior,
                        data = data,
                        sampler = sampler,
                        n_iter = n_iter,
                        n_thin = 1,
                        save = TRUE, prog = 0,
                        upper_tau2 = upper_tau2)
  })

  draws <- sapply(mcmc_out$sims, function(x) x[[param]])

  if (isTRUE(sampler[[param]]$logscale)) {
    draws_ess <- log(draws) |> coda::as.mcmc()
  } else {
    draws_ess <- draws |> coda::as.mcmc()
  }

  ESS <- min(n_iter, coda::effectiveSize(draws_ess))
  out_df <- data.frame(nEval = sapply(mcmc_out$extras, function(x) {
                                          x$nEvaluations
                                        }) |> sum(),
                    EffSamp = ESS,
                    userTime = time_out['user.self'],
                    sysTime = time_out['sys.self'],
                    elapsedTime = time_out['elapsed'])

  list(timing = out_df, draws = draws, state = mcmc_out$state)
}
