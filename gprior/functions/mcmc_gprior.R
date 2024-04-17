# source("full_conditionals.R")

mcmc_gprior <- function(state, prior, data, sampler, n_iter, save = TRUE, prog = 0) {

  ## state is a list with: beta, psi, g, and latent_s (latent slice for g), and iter
  ## prior is a list with: beta_0, a_0, b_0, g_max
  ## data is a list with: X, y, p = ncol(X), n = length(y) = nrow(X),
  ##    XtX, inv_chol_XtX (where Chol is upper triangular), logdet_XtX, beta_mle

  if (save) {
    sims <- rep(list(state), n_iter)
    extras <- rep(list(), n_iter)
  } else {
    sims <- NULL
    extras <- NULL
  }

  for (i in 1:n_iter) {

    state$beta <- update_beta(state = state, prior = prior, data = data)

    state$psi <- update_psi(state = state, prior = prior, data = data)

    tmp <- update_g(state = state, prior = prior, data = data,
                    sampler = sampler[["g"]])
    state <- tmp$state

    state$iter <- state$iter + 1

    if (save) {
      sims[[i]] <- state
      extras[[i]] <- tmp[["extras"]]
    }

    if (prog > 0) {
      if (i %% prog == 0) {
        cat("Iter ", i, " of ", n_iter, "\n")
        timestamp()
      }
    }

  }

  list(state = state, prior = prior, data = data, sampler = sampler,
       save = save, sims = sims, extras = extras, n_iter = n_iter)
}


time_gprior <- function(state, prior, data, sampler, n_iter) {

  require("coda")

  time_out <- system.time({
    mcmc_out <- mcmc_gprior(state = state,
                            prior = prior,
                            data = data,
                            sampler = sampler,
                            n_iter = n_iter,
                            save = TRUE, prog = 0)
  })

  draws_g <- sapply(mcmc_out$sims, function(x) x$g)

  if (sampler[["g"]]$logG) {
    draws_ess <- log(draws_g) |> coda::as.mcmc()
  } else {
    draws_ess <- draws_g |> coda::as.mcmc()
  }

  ESS <- min(n_iter, coda::effectiveSize(draws_ess))
  out_df <- data.frame(nEval = sapply(mcmc_out$extras, function(x) {
                                          x$nEvaluations
                                        }) |> sum(),
                    EffSamp = ESS,
                    userTime = time_out['user.self'],
                    sysTime = time_out['sys.self'],
                    elapsedTime = time_out['elapsed'])

  list(timing = out_df, draws = draws_g, state = mcmc_out$state)
}
