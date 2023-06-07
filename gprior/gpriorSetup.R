## gprior setup function
## author: Sam Johnson

defaultW <- getOption("warn") 

options(warn = -1) 

library(mvtnorm)
library(cucumber)
library(doParallel)
library(mcmc)
library(coda)

source('pseudoCauchyFunctions.R')
source('formatingFunctions.R')

Nburnin <- 1000
Nthin <- 2
Nchains <- 5
Nsamples <- 50000

print(paste0('Total Samples after burnin and thinning is: ', ((Nsamples - Nburnin)/Nthin) * Nchains))

data(mtcars)

attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)

# beta ~ mvn(beta_0, g / psi * inv(XtX)), using the covariance parametrization
yty = t(y) %*% y
Xt <- t(X)
inv_XtX <- solve(Xt %*% X)
beta_mle <- solve(Xt %*% X, Xt %*% y)
beta_0 <- rep(0, ncol(X))
beta <- beta_0

# psi ~ gamma(a, b)
a_0 <- 5
b_0 <- 1
psi <- a_0 / b_0

# g ~ uniform(0, g_max)
g_max <- 2 * ncol(X)^2
g <- g_max / 2


# pre allocating memory
samples <- list(beta = matrix(0.0, nrow = Nsamples, ncol = length(beta)), psi = numeric(Nsamples),
                g = numeric(Nsamples), u = numeric(Nsamples), nEval = numeric(Nsamples),
                time = numeric(1))
# setting up chain storage
chainSamples <- vector('list', length = Nchains)
chainSamples <- lapply(chainSamples, \(list) samples)


## functions

# function needed in mcmc

log_target <- function(g) {
  sapply(g, FUN = \(g) {
  beta_cov <- g / psi * inv_XtX
  dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)}
  )
}

target <- function(g) {
  sapply(g, FUN = \(g) {
  beta_cov <- g / psi * inv_XtX
  dmvnorm(beta, beta_0, beta_cov) * ifelse(0 < g && g < g_max, 1, 0)}
  )
}

update_parameters <- function() {
  q <- g / (1 + g)
  beta_cov <<- q / psi * inv_XtX
  beta_mean <<- q * beta_mle + (1 - q) * beta_0
  beta <<- rmvnorm(1, beta_mean, beta_cov)
  a_n <<- a_0 + n / 2
  residuals <- y - X %*% t(beta)
  b_n <<- b_0 + 0.5 * sum(residuals^2)
  psi <<- rgamma(1, a_n, b_n)
}

saving_updates <- function(method = 'Other') {
  chainSamples[[chain]]$beta[i,] <<- beta
  chainSamples[[chain]]$psi[i] <<- psi
  chainSamples[[chain]]$g[i] <<- g
  chainSamples[[chain]]$nEval[i] <<- nEval
  if(method == 'Transform') chainSamples[[chain]]$u[i] <<- u
}

burnin_thinning <- function(burnin = Nburnin, thin = Nthin, method = 'Other') {
  if(burnin < 1) stop('burnin has to be greater than or equal to 1')
  chainSamples[[chain]]$beta <<- chainSamples[[chain]]$beta[-c(1:burnin),]
  chainSamples[[chain]]$beta <<- chainSamples[[chain]]$beta[LaplacesDemon::Thin(1:nrow(chainSamples[[chain]]$beta), By = thin),]
  chainSamples[[chain]]$psi <<- LaplacesDemon::Thin(chainSamples[[chain]]$psi[-c(1:burnin)], By = thin)
  chainSamples[[chain]]$g <<- LaplacesDemon::Thin(chainSamples[[chain]]$g[-c(1:burnin)], By = thin)
  chainSamples[[chain]]$nEval <<- LaplacesDemon::Thin(chainSamples[[chain]]$nEval[-c(1:burnin)], By = thin)
  if(method == 'Transform') chainSamples[[chain]]$u <<- LaplacesDemon::Thin(chainSamples[[chain]]$u[-c(1:burnin)], By = thin)
}

save_time <- function() {
  chainSamples[[chain]]$time <<- time['user.self']
}

# function dependent on sampler

# rand walk function
# proposed draw
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

## transform sampler

# get samples using stepping out procedure

approx_samples <- function(x, lf, samples, w) {
  draws <- numeric(samples)
  draws[1] <- x
  for( i in 2:samples ) {
    if(sum(is.na(draws)) != 0) browser()
    draws[i] <- cucumber::slice_sampler_stepping_out(draws[i-1], target = lf, w = w, log = TRUE, max = Inf)$x
  }
  draws
}


# setting up parallel processing
cl <- makeCluster(Nchains)
registerDoParallel(cl)

options(warn = defaultW)