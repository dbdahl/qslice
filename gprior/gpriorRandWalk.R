## g prior sampling using rand walk
## Author: Sam Johnson

library(mvtnorm)
library(cucumber)
library(doParallel)
library(mcmc)
library(coda)

data(mtcars)

attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)

nChains <- 5
cl <- makeCluster(nChains)
registerDoParallel(cl)

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
      # n.accept <- n.accept+1
    }
  }
  list(x = int.x, accept = accept)
}


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
samples <- list(beta = matrix(0.0, nrow = n_samples, ncol = length(beta)), psi = numeric(n_samples), g = numeric(n_samples), accept = numeric(n_samples), time = numeric(1))
# setting up chain storage
chainSamples <- vector('list', length = nChains)
chainSamples <- lapply(chainSamples, \(list) samples)

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
      temp <- randWalk(int.x = g, lf = \(g) {
        beta_cov <- g / psi * inv_XtX
        dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
      }, c = 30)
      g <- temp$x
      chainSamples[[chain]]$beta[i,] <- beta
      chainSamples[[chain]]$psi[i] <- psi
      chainSamples[[chain]]$g[i] <- g
      chainSamples[[chain]]$accept[i] <- temp$accept
    }
  })
  chainSamples[[chain]]$time <- time['user.self']
}

# acceptance rate
sapply(chainSamples, \(list) mean(list$accept))


# evaluation of samples

traceplot <- function(matrix, mainTitle = NULL) {
  n <- ncol(matrix)
  # finding eff sample size of samples
  effSize <- round(sum(apply(matrix, 2, effectiveSize)))
  # finding Rhat
  chainList <- mcmc.list()
  for (i in 1:n) chainList[[i]] <- mcmc(matrix[,i])
  Rhat <- round(gelman.diag(chainList)$psrf[1,1],3)
  # making trace plot
  color <- RColorBrewer::brewer.pal(n, 'Set2')
  plot(matrix[,1], type = 'l', col = color[i], main = paste0(mainTitle, ' ESS=',effSize, ' Rhat=',Rhat),
       ylab = '')
  for( i in 2:n ) {
    lines(matrix[,i], type = 'l', col = color[i])
  }
}


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
sapply(chainSamples, \(list) list$time)


saveRDS(chainSamples, file = 'data/randWalk.rds')
