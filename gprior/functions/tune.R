tune <- function(state, prior, data, sampler, bnds_init,
                 n_iter = 1000,
                 n_grid = 5,
                 n_rep = 3,
                 n_rounds = 3,
                 range_frac = 0.5,
                 verbose = TRUE) {

  require("coda")
  stopifnot(sampler$g$type %in% c("rw", "stepping", "latent"))

  esps <- matrix(NA, nrow = n_grid, ncol = n_rep)
  hstry <- list()
  bnds_now <- bnds_init
  range_now <- diff(bnds_now)

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

        if (sampler$g$type == "rw") {
          sampler$g$c = vals[i]
        } else if (sampler$g$type == "stepping") {
          sampler$g$w = vals[i]
        } else if (sampler$g$type == "latent") {
          sampler$g$rate = vals[i]
        }

        ## collect timing info
        mc_time <- time_gprior(state = state, prior = prior, data = data,
                               sampler = sampler, n_iter = n_iter)
        state <- mc_time$state

        esps[i,j] <- mc_time$timing$EffSamp / mc_time$timing$userTime

      }
    }

    ## find max
    esps_means <- rowMeans(esps)
    indx_max <- which.max(esps_means)

    ## collect history
    hstry[[rr]] <- list(vals = vals, esps = esps, means = esps_means,
                        opt = vals[indx_max])

    if (isTRUE(verbose)) {
      cat("Round", rr, "\n", "values/ESPS:",
          paste(vals, round(esps_means), collapse = "; "),
          "\n", "Selection:",
          vals[indx_max], "at", round(esps_means[indx_max]),
          "\n")
    }

    ## set up subsequent grid
    if (rr <= n_rounds) {
      bnds_now[1] <- vals[indx_max] - range_now * range_frac / 2.0
      bnds_now[2] <- vals[indx_max] + range_now * range_frac / 2.0
      range_now <- diff(bnds_now)
    }

  }

  ## final selection
  val_opt <- vals[indx_max]

  if (isTRUE(verbose)) {
    cat("Final selection:",
        val_opt, "at", round(esps_means[indx_max]),
        "\n")
  }

  if (sampler$g$type == "rw") {
    sampler$g$c = val_opt
  } else if (sampler$g$type == "stepping") {
    sampler$g$w = val_opt
  } else if (sampler$g$type == "latent") {
    sampler$g$rate = val_opt
  }

  list(val_opt = val_opt, hstry = hstry, state = state,
       sampler = sampler, n_iter = n_iter)

}
