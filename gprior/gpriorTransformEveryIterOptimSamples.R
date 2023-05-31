## g prior sampling using transform slice psuedo updated at every iteration using optim based on samples
## Author: Sam Johnson

library(mvtnorm)
library(cucumber)
library(doParallel)
library(mcmc)
library(coda)

source('pseudoCauchyFunctions.R')
source('formatingFunctions.R')

data(mtcars)

attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)

nChains <- 5

# fm <- lm(y ~ -1 + X)
# summary(fm)

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

# number of samples to be taken
n_samples <- 10000
# pre allocating memory
samples <- list(beta = matrix(0.0, nrow = n_samples, ncol = length(beta)), psi = numeric(n_samples), g = numeric(n_samples), u = numeric(n_samples), time = numeric(1))
# setting up chain storage
chainSamples <- vector('list', length = nChains)
chainSamples <- lapply(chainSamples, \(list) samples)

# fitting a pseudo target
beta <- beta + c(0.04, 0.46, -0.38, 0.09, -0.36, 0.09, 0.10, 0.20, 0.41, -0.27)
f <- \(g) {
  beta_cov <- g / psi * inv_XtX
  exp(dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf))
}

xx <- seq(from = 0.01, to = 200, length.out = 1000)
yy <- sapply(xx,f)

# plot(xx,yy)

# fitting a pseudo target
approxSamples <- numeric(5000)

for (i in seq_along(approxSamples)) {
  q <- g / (1 + g)
  beta_cov <- q / psi * inv_XtX
  beta_mean <- q * beta_mle + (1 - q) * beta_0
  beta <- rmvnorm(1, beta_mean, beta_cov)
  a_n <- a_0 + n / 2
  residuals <- y - X %*% t(beta)
  b_n <- b_0 + 0.5 * sum(residuals^2)
  psi <- rgamma(1, a_n, b_n)
  g <- cucumber::slice_sampler_stepping_out(g, \(g) {
    beta_cov <- g / psi * inv_XtX
    dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
  }, w = 50, log = TRUE)$x
  
  approxSamples[i] <- g
}

psuedoFit <- opt_Cauchy_auc_data(approxSamples, lb = 0)
psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)

cl <- makeCluster(nChains)
registerDoParallel(cl)

# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
  time <- system.time({
    for (i in seq_len(n_samples)) {
      q <- g / (1 + g)
      beta_cov <- q / psi * inv_XtX
      beta_mean <- q * beta_mle + (1 - q) * beta_0
      beta <- rmvnorm(1, beta_mean, beta_cov)
      a_n <- a_0 + n / 2
      residuals <- y - X %*% t(beta)
      b_n <- b_0 + 0.5 * sum(residuals^2)
      psi <- rgamma(1, a_n, b_n)
      ##
      if(i %% 1 == 1000) {
        newg <- chainSamples[[chain]]$g[1:i]
        if(min(newg) < 0) browser(text = 'A g is less than 0')
        psuedoFit <- opt_Cauchy_auc_data(newg, lb = 0)
        psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = psuedoFit$lb, ub = psuedoFit$ub)
      }
      ##
      temp <- cucumber::slice_sampler_transform(x = g, target = \(g) {
        beta_cov <- g / psi * inv_XtX
        dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
      }, pseudo_log_pdf = psuedoTarget$pseudo_log_pdf, pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)
      ##
      # f <- \(g) {
      #   beta_cov <- g / psi * inv_XtX
      #   exp(dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf))
      # }
      # psuedoFit <- lapproxt(f = f, init = 10, lb = 0)
      # psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = 0, name = 'Laprox')
      ##
      g <- temp$x
      chainSamples[[chain]]$beta[i,] <- beta
      chainSamples[[chain]]$psi[i] <- psi
      chainSamples[[chain]]$g[i] <- g
      chainSamples[[chain]]$u[i] <- temp$u
    }
    chainSamples[[chain]]$beta <- chainSamples[[chain]]$beta[-c(1:Nburnin),]
    chainSamples[[chain]]$beta <- chainSamples[[chain]]$beta[LaplacesDemon::Thin(1:nrow(chainSamples[[chain]]$beta), By = Nthin),]
    chainSamples[[chain]]$psi <- LaplacesDemon::Thin(chainSamples[[chain]]$psi[-c(1:Nburnin)], By = Nthin)
    chainSamples[[chain]]$g <- LaplacesDemon::Thin(chainSamples[[chain]]$g[-c(1:Nburnin)], By = Nthin)
    chainSamples[[chain]]$u <- LaplacesDemon::Thin(chainSamples[[chain]]$u[-c(1:Nburnin)], By = Nthin)
  })
  chainSamples[[chain]]$time <- time['user.self']
}

# evaluation of samples

saveRDS(chainSamples, file = 'data/transformEveryIterOptimSamples.rds')


print('finished')

