
############## Random Walk ################

random_walk_sampler <- function(n_iter, lf, support, x_0, c) {
  draws <- numeric(n_iter + 1)
  draws[1] <- x_0
  int.x <- x_0
  n.accept <- 0
  for(i in 2:(n_iter + 1)){
    ## proposed draw
    x.dot <- rnorm(1, mean = int.x, sd = c)
    if (x.dot >= support[1] && x.dot <= support[2]) {

      logr <- lf(x.dot) - lf(int.x)

      u <- runif(1, 0, 1)

      if(log(u) < logr){
        int.x <- x.dot
        n.accept <- n.accept + 1
      }
    }
    draws[i] <- int.x
  }
  list(draws = draws[-1], counter = 2*n_iter, n.accept = n.accept, n_samp = n_iter)
}


############ Stepping Out Eval #############


# function to evaluate stepping out procedure

stepping_out_sampler <- function(n_iter, lf, x_0, w, max) {
  counter <- 0
  draws <- numeric(n_iter + 1)
  draws[1] <- x_0

  for ( i in 2:(n_iter + 1)) {
    out <- slice_stepping_out(x = draws[i-1], log_target = lf, w = w, max = max)
    draws[i] <- out$x
    counter <- counter + out$nEvaluations
  }

  list(draws = draws[-1], counter = counter, n_samp = n_iter)
}


################ GESS ##############

## generalized elliptical slice sampler (Nishihara 2014)

gess_sampler <- function(n_iter, lf, x_0, mu, sigma, degf) {
  counter <- 0
  draws <- numeric(n_iter + 1)
  draws[1] <- x_0

  for ( i in 2:(n_iter + 1) ) {
    out <- slice_genelliptical(x = draws[i-1], log_target = lf,
                               mu = mu, sigma = sigma, df = degf)

    draws[i] <- out$x
    counter <- counter + out$nEvaluations
  }

  list(draws = draws[-1], counter = counter, n_samp = n_iter)
}


############ Latent Eval ############

# creating a function to evaluate the latent slice sampler

latent_sampler <- function(n_iter, lf, x_0, s_0, rate) {
  counter <- 0
  draws <- latents <- numeric(n_iter + 1)
  draws[1] <- x_0
  latents[1] <- s_0

  for ( i in 2:(n_iter + 1) ) {
    out <- slice_latent(x = draws[i-1], s = latents[i-1],
                        log_target = lf, rate = rate)
    draws[i] <- out$x
    latents[i] <- out$s
    counter <- counter + out$nEvaluations
  }

  list(draws = draws[-1], counter = counter, n_samp = n_iter)
}



############## Quantile Slice Eval ################

quantile_sampler <- function(n_iter, lf, x_0, pseudo_lpdf, pseudo_inv_cdf) {

  counter <- 0
  draws <- numeric(n_iter + 1)
  Tdraws <- numeric(n_iter + 1)
  draws[1] <- x_0
  Tdraws[1] <- 0

  for ( i in 2:(n_iter + 1) ) {
    out <- slice_quantile(x = draws[i-1],
                          log_target = lf, pseudo_log_pdf = pseudo_lpdf,
                          pseudo_inv_cdf = pseudo_inv_cdf)
    draws[i] <- out$x
    Tdraws[i] <- out$u
    counter <- counter + out$nEvaluations
  }

  list(draws = draws[-1], Tdraws = Tdraws[-1], counter = counter, n_samp = n_iter)
}


############## Independence Metropolis Hastings ################

IMH_sampler <- function(n_iter, lf, x_0, pseudo) {
  draws <- numeric(n_iter + 1)
  draws[1] <- x_0
  n.accept <- 0
  for(i in 2:(n_iter + 1)){
    tmp <- imh_pseudo(x = draws[i-1], log_target = lf, pseudo = pseudo)
    draws[i] <- tmp$x
    n.accept <- tmp$accpt
  }
  list(draws = draws[-1], counter = 2*n_iter, n.accept = n.accept, n_samp = n_iter)
}



### universal timer

# function to evaluate Random Walk
sampler_time_eval <- function(type,
                              n_iter,
                              lf_func,
                              support,
                              x_0,
                              settings) {

  if (type == "rw") {

    cc <- settings[1, "c"]

    time <- system.time({
      mcmc_out <- random_walk_sampler(n_iter = n_iter,
                                      lf = lf_func,
                                      support = support,
                                      x_0 = x_0,
                                      c = cc)
    })

  } else if (type == "stepping") {

    w <- settings[1, "w"]

    time <- system.time({
      mcmc_out <- stepping_out_sampler(n_iter = n_iter,
                                       lf = lf_func,
                                       x_0 = x_0,
                                       w = w,
                                       max = Inf)
    })

  } else if (type == "gess") {

    loc <- settings[1, "loc"]
    sc <- settings[1, "sc"]
    degf <- settings[1, "degf"]

    time <- system.time({
      mcmc_out <- gess_sampler(n_iter = n_iter,
                               lf = lf_func,
                               x_0 = x_0,
                               mu = loc,
                               sigma = sc,
                               degf = degf)
    })

  } else if (type == "latent") {

    s_0 <- settings[1, "s_init"]
    rate <- settings[1, "rate"]

    time <- system.time({
      mcmc_out <- latent_sampler(n_iter = n_iter,
                                 lf = lf_func,
                                 x_0 = x_0,
                                 s_0 = s_0,
                                 rate = rate)
    })

  } else if (type == "Qslice") {

    time <- system.time({
      mcmc_out <- quantile_sampler(n_iter = n_iter,
                                   lf = lf_func,
                                   x_0 = x_0,
                                   pseudo_lpdf = settings$ld,
                                   pseudo_inv_cdf = settings$q)
    })

  } else if (type == "imh") {

    time <- system.time({
      mcmc_out <- IMH_sampler(n_iter = n_iter,
                              lf = lf_func,
                              x_0 = x_0,
                              pseudo = pseudo)
    })

  }

  ESS <- min(n_iter, coda::effectiveSize(coda::as.mcmc(mcmc_out$draws)))
  temp_tbl <- data.frame(nEval = mcmc_out$counter, EffSamp = ESS,
                         userTime = time['user.self'],
                         sysTime = time['sys.self'],
                         elapsedTime = time['elapsed'])

  temp_tbl$Draws <- list(mcmc_out$draws)
  temp_tbl
}
