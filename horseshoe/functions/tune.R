tune <- function(state, prior, data, sampler, param, bnds_init,
                 n_iter = 1000,
                 n_grid = 3,
                 n_rep = 2,
                 n_rounds = 5,
                 range_frac = 0.67,
                 verbose = TRUE) {

  require("coda")
  stopifnot(sampler[[param]]$type %in% c("rw", "stepping", "latent"))

  esps <- matrix(NA, nrow = n_grid, ncol = n_rep)
  hstry <- list()
  bnds_now <- bnds_init
  range_now <- diff(bnds_now)

  esps_means_running <- numeric(0)
  vals_running <- numeric(0)
  lxx <- numeric(0)
  lyy <- numeric(0)

  for (rr in 1:n_rounds) {

    vals <- seq(bnds_now[1], bnds_now[2], length = n_grid)

    if (any(vals <= 0.0)) {
      indx_neg <- which(vals <= 0.0)
      min_pos <- vals[-indx_neg][1]
      replacement_vals <- seq(0, min_pos, length = length(indx_neg) + 2)
      vals[indx_neg] <- replacement_vals[2:(length(replacement_vals) - 1)]
      bnds_now[1] <- vals[1]
      range_now <- diff(bnds_now)
    }

    for (i in 1:length(vals)) {
      for (j in 1:n_rep) {

        if (sampler[[param]]$type == "rw") {
          sampler[[param]]$c = vals[i]
        } else if (sampler[[param]]$type == "stepping") {
          sampler[[param]]$w = vals[i]
        } else if (sampler[[param]]$type == "latent") {
          sampler[[param]]$rate = vals[i]
        }

        ## collect timing info
        mc_time <- time_hs(state = state, prior = prior, data = data,
                           sampler = sampler, n_iter = n_iter, param = param)
        state <- mc_time$state

        esps[i,j] <- mc_time$timing$EffSamp / mc_time$timing$userTime

      }
    }

    ## collect running results
    esps_means <- rowMeans(esps)
    esps_means_running <- c(esps_means_running, esps_means)
    vals_running <- c(vals_running, vals)

    ## find max

    # first, try a quadratic regression on all results so far

    lxx <- c(lxx, rep(log(vals), each = n_rep))
    lyy <- c(lyy, c(t(log(esps))))

    mod <- lm(lyy ~ lxx + I(lxx^2))
    betas <- coef(mod)

    if (betas[3] < 0.0) { # proceed with quadratic

      lx_max <- -0.5 * betas[2] / betas[3]
      opt <- exp(lx_max)

    } else {

      opt <- vals_running[which.max(esps_means_running)]

    }

    ## collect history
    hstry[[rr]] <- list(vals = vals, esps = esps, means = esps_means,
                        opt = opt, lvals = lxx, lesps = lyy)

    if (isTRUE(verbose)) {
      cat("Round", rr, "\n", "values/ESPS:",
          paste(vals, round(esps_means, 2), collapse = "; "),
          "\n", "Selection:",
          opt,
          "\n")
    }

    ## set up subsequent grid
    if (rr <= n_rounds) {
      bnds_now[1] <- opt - range_now * range_frac / 2.0
      bnds_now[2] <- opt + range_now * range_frac / 2.0
      range_now <- diff(bnds_now)
    }

  }

  ## final selection

  if ( (opt < min(vals_running)) || (opt > max(vals_running)) ) { # if solution out of bounds, don't extrapolate
    val_opt <- vals_running[which.max(esps_means_running)]
  } else {
    val_opt <- opt
  }

  if (isTRUE(verbose)) {
    cat("Final selection:",
        val_opt,
        "\n")
  }

  if (sampler[[param]]$type == "rw") {
    sampler[[param]]$c = val_opt
  } else if (sampler[[param]]$type == "stepping") {
    sampler[[param]]$w = val_opt
  } else if (sampler[[param]]$type == "latent") {
    sampler[[param]]$rate = val_opt
  }

  list(val_opt = val_opt, lvals = lxx, lesps = lyy,
       hstry = hstry, state = state,
       sampler = sampler, n_iter = n_iter)
}
