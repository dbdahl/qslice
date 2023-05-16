## g prior sampling using transform slice psuedo updated at every interation using laplace
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
cl <- makeCluster(nChains)
registerDoParallel(cl)

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
samples <- list(beta = matrix(0.0, nrow = n_samples, ncol = length(beta)), psi = numeric(n_samples), g = numeric(n_samples), time = numeric(1))
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

plot(xx,yy)

# this doesn't work because there is no maximum
psuedoFit <- lapproxt(f, 10, lb = 0)
psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = 0, name = 'Laprox')

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
      f <- \(g) {
        beta_cov <- g / psi * inv_XtX
        exp(dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf))
      }
      psuedoFit <- lapproxt(f, 10, lb = 0)
      psuedoTarget <- pseudo_Cauchy(loc = psuedoFit$loc, sc = psuedoFit$sc, lb = 0, name = 'Laprox')
      browser()
      g <- cucumber::slice_sampler_transform(x = g, target = \(g) {
        beta_cov <- g / psi * inv_XtX
        dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
      }, pseudo_log_pdf = psuedoTarget$pseudo_log_pdf, pseudo_inv_cdf = psuedoTarget$pseudo_inv_cdf)$x
      chainSamples[[chain]]$beta[i,] <- beta
      chainSamples[[chain]]$psi[i] <- psi
      chainSamples[[chain]]$g[i] <- g
    }
  })
  chainSamples[[chain]]$time <- time['user.self']
}

# evaluation of samples

# g parameter
gSamples <- sapply(chainSamples, \(list) list$g)

traceplot(gSamples,'g')
plot(density(c(gSamples)), main = 'g')

# psi parameter
psiSamples <- sapply(chainSamples, \(list) list$psi)

traceplot(psiSamples, 'psi')
plot(density(c(psiSamples)), main = 'psi')

# beta parameters
betaSamples <- lapply(1:ncol(X), \(i) sapply(chainSamples, \(list) list$beta[,i]))

par(mfrow = c(3,4))
sapply(betaSamples, traceplot)

par(mfrow = c(3,4))
sapply(betaSamples, \(list) plot(density(c(list)), main = 'beta'))

# time
time <- sapply(chainSamples, \(list) list$time)


saveRDS(chainSamples, file = 'data/transformEveryIterLaplace.rds')




