# functions use to update the gprior
# author: Sam Johnson


# rand walk function
randWalk <- function(int.x, lf, c, support = c(-Inf, Inf)) {
  x.dot <- rnorm(1, int.x, c)
  if(x.dot >= support[1] && x.dot <= support[2]){
    
    logr <- lf(x.dot) - lf(int.x)
    accept <- 0
    u <- runif(1, 0, 1)
    if(log(u) < logr){
      int.x <- x.dot
      accept <- 1
    }
  }
  list(x = int.x, accept = accept)
}

# independence mh
independence <- function(int.x, psuedoTarget, lf) {
  
  x.dot <- psuedoTarget$q(runif(1))
  logr <- (lf(x.dot) + psuedoTarget$ld(int.x)) - (lf(int.x) + psuedoTarget$ld(x.dot))
  accept <- 0
  u <- runif(1, 0, 1)
  if(log(u) < logr){
    int.x <- x.dot
    accept <- 1
  }
  list(x = int.x, accept = accept)
}


update_g <- function(g, psi, beta, list, i) {
  updateGEnv <- environment()
  list2env(list, envir = updateGEnv)
  # defining the target and log target inside this function
  # log target
  log_target <- function(g) {
    sapply(g, FUN = \(g) {
      beta_cov <- g / psi * inv_XtX
      dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
    }
    )
  }
  # target
  target <- function(g) {
    sapply(g, FUN = \(g) {
      beta_cov <- g / psi * inv_XtX
      dmvnorm(beta, beta_0, beta_cov) * ifelse(0 < g && g < g_max, 1, 0)
    }
    )
  }
  ##### h
  h <- function(g) {
    exp(log_target(pseudoTarget$pseu$q(g)) - pseudoTarget$pseu$ld(pseudoTarget$pseu$q(g)))
  }
  ####
  temp <- if(method == 'SteppingOut') {
    temp <- cucumber::slice_sampler_stepping_out(x = g, target = log_target,
                                                 w = w, log = TRUE, max = Inf)
    list(g = temp$x, u = 0, nEval = temp$nEvaluations)
  } else if (method == 'Latent') {
    temp <- cucumber::slice_sampler_latent(x = g, s = s, target = log_target,
                                           rate = rate, log = TRUE)
    list(g = temp$x, u = 0, nEval = temp$nEvaluations)
  } else if (method == 'GESS') {
    temp <- cucumber::slice_sampler_generalized_elliptical(x = g, target = log_target,
                                                           mu = mu, sigma = sigma, df = df, log = TRUE)
    list(g = temp$x, u = 0, nEval = temp$nEvaluations)
  } else if (method == 'RandWalk') {
    temp <- randWalk(int.x = g, lf = log_target, c = c)
    list(g = temp$x, u = 0, nEval = temp$accept)
  } else if (method == 'Independence') {
    temp <- independence(int.x = g, psuedoTarget = pseudoTarget$pseu, lf = log_target)
    list(g = temp$x, u = 0, nEval = temp$accept)
  } else if (method == 'Transform') {
    temp <- cucumber::slice_sampler_transform(x = g, target = log_target,
                                              pseudo_log_pdf = pseudoTarget$pseu$ld,
                                              pseudo_inv_cdf = pseudoTarget$pseu$q)
    list(g = c(temp$x), u = temp$u, nEval = temp$nEvaluations)
  }
  temp
}

# create pseudo target
create_pseudo <- function(psi, beta, pseudoList) {
  createPseudoEnv <- environment()
  list2env(pseudoList, envir = createPseudoEnv)
  # checks
  if(!(pseudoType %in% c('Laplace','Auto','OptimSamples','Optim','OptimSamplesAUC','OptimAUC'))) {
    stop('pseudoType must be one of these: Laplace, Auto, OptimSamples, Optim, OptimSamplesAUC, OptimAUC')
  }
  # creating the functions that the pseudo target will be make from
  # log target
  log_target <- function(g) {
    sapply(g, FUN = \(g) {
      beta_cov <- g / psi * inv_XtX
      dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
    }
    )
  }
  # target
  target <- function(g) {
    sapply(g, FUN = \(g) {
      beta_cov <- g / psi * inv_XtX
      dmvnorm(beta, beta_0, beta_cov) * ifelse(0 < g && g < g_max, 1, 0)
    }
    )
  }
  # fitting the pseudo target
  tempPseudo <- if (pseudoType == 'Laplace') {
    c <- 0.5 * psi * (beta-beta_0)%*%XtX%*%t(beta-beta_0)
    loc <- (2*c)/(n)
    sigma <- function(g) {
      x1 <- (2*c)/(g^3)
      x2 <- (n)/(2*g^2)
      x1x2 <- x1 - x2
      sqrt(1/x1x2)
    }
    pseu <- pseudo_t_list(loc = loc, sc = sigma(loc), degf = 1, lb = 0)#lapproxt(target, 10, lb = 0, maxub = 200)
    pseudoTarget <- list(pseu = pseu, name = 'laplace')
    pseudoTarget
  } else if (pseudoType == 'Auto') {
    pseudoTarget <- fit_trunc_Cauchy(approxSamples, lb = 0)
    pseudoTarget
  } else if (pseudoType == 'OptimSamples') {
    pseudoTarget <- opt_t(samples = approxSamples, lb = 0, nbins = 30, type = optimType, coeffs = c(1,0))
    pseudoTarget
  } else if (pseudoType == 'OptimSamplesAUC') {
    pseudoTarget <- opt_t(samples = approxSamples, lb = 0, nbins = 30, type = optimType, use_meanSliceWidth = FALSE, coeffs = c(1,0))
    pseudoTarget
  } else if (pseudoType == 'Optim') {
    truth = list(ld = log_target,
                 t = "Full Conditional of g",
                 lb = 0,
                 ub = 200)
    truth$dld = function(x) numDeriv::grad(truth$ld, x=x)
    pseudoTarget <- opt_t(target = truth, lb = 0, nbins = 500, type = optimType, coeffs = c(1,0))
    pseudoTarget
  } else if (pseudoType == 'OptimAUC') {
    truth = list(ld = log_target,
                 t = "Full Conditional of g",
                 lb = 0,
                 ub = 200)
    truth$dld = function(x) numDeriv::grad(truth$ld, x=x)
    pseudoTarget <- opt_t(target = truth, lb = 0, nbins = 500, type = optimType, use_meanSliceWidth = FALSE, coeffs = c(1,0))
    pseudoTarget
  }
  tempPseudo
}

# fixedBeta <- c(-0.02742211, 0.25728077, -0.23017904, 0.06894651, -0.56840819,
#                0.23139456, 0.02402770,  0.19669150, 0.07542407, -0.05148203)
# fixedPsi <- 5.264618

gprior_sampler <- function(list) {
  # browser()
  # getting the environment of the function
  env <- environment()
  list2env(list, env)
  if(!(method %in% c('SteppingOut','Latent','GESS','RandWalk','Transform','Independence'))) {
    stop('Method must be one of the following: SteppingOut, Latent, GESS, RandWalk, Transform','Independence')
  }
  # creating a list to store all the variables
  sampleList <- list(beta = matrix(0.0, nrow = Nsamples, ncol = length(beta)), psi = numeric(Nsamples),
                     g = numeric(Nsamples), u = numeric(Nsamples), nEval = numeric(Nsamples),
                     time = numeric(1))
  # setting initial values
  time <- system.time({
    for(i in 2:Nsamples) {
      # updating beta and psi
      q <- g / (1 + g)
      beta_cov <- q / psi * inv_XtX
      beta_mean <- q * beta_mle + (1 - q) * beta_0
      beta <- rmvnorm(1, beta_mean, beta_cov) # t(fixedBeta)
      a_n <- a_0 + n #a_0 + n / 2
      residuals <- y - X %*% t(beta)
      betaDiff <- beta - beta_0
      b_n <- 0.5 * sum(residuals^2) + 0.5 * betaDiff %*% (XtX/g) %*% t(betaDiff) + b_0 #b_0 + 0.5 * sum(residuals^2)
      psi <- rgamma(1, a_n, b_n) #fixedPsi 
      # creating a new pseudo target
      if(method == 'Transform' || method == 'Independence') {
        if (i == 2) {
          list$pseudoTarget <- create_pseudo(psi = psi, beta = beta, pseudoList = list)
        }
        if(is.numeric(everyIter) && i != 2) {
          if( i %% everyIter == 0) {
            list$pseudoTarget <- create_pseudo(psi = psi, beta = beta, pseudoList = list)
          }
        }
      }
      # update g
      tempList <- update_g(g = g, psi = psi, beta = beta, list = list, i = i)
      list2env(tempList, envir = env)
      # saving the updates
      sampleList$beta[i,] <- beta
      sampleList$psi[i] <- psi
      sampleList$g[i] <- g
      sampleList$u[i] <- u
      sampleList$nEval[i] <- nEval
    }
  })
  # thinning and burnin 
  sampleList$beta <- sampleList$beta[-c(1:Nburnin),]
  sampleList$beta <- sampleList$beta[LaplacesDemon::Thin(1:nrow(sampleList$beta), By = Nthin),]
  sampleList$psi <- sampleList$psi[-c(1:Nburnin)] |> LaplacesDemon::Thin(By = Nthin)
  sampleList$g <- sampleList$g[-c(1:Nburnin)] |> LaplacesDemon::Thin(By = Nthin)
  sampleList$nEval <- sampleList$nEval[-c(1:Nburnin)] |> LaplacesDemon::Thin(By = Nthin)
  sampleList$u <- sampleList$u[-c(1:Nburnin)] |> LaplacesDemon::Thin(By = Nthin)
  sampleList$time <- time
  
  sampleList
}
