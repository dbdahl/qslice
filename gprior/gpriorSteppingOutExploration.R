## exploring stepping out for the g prior
## Author: Sam Johnson

library(mvtnorm)
library(cucumber)
library(doParallel)
library(mcmc)
library(coda)

source('formatingFunctions.R') # needed to create the traceplots

slice_sampler_stepping_out <- function (x, target, w, max = 0, log = TRUE) 
{
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
  }
  fx <- f(x)
  y <- if (isTRUE(log)) {
    log(runif(1)) + fx
  }
  else {
    runif(1) * fx
  }
  L <- x - runif(1) * w
  R <- L + w
  if (!is.finite(max)) {
    while (y < f(L)) L <- L - w
    while (y < f(R)) R <- R + w
  }
  else if (max > 0) {
    J <- floor(runif(1) * max)
    K <- max - 1 - J
    while (J > 0 && y < f(L)) {
      L <- L - w
      J <- J - 1
    }
    while (K > 0 && y < f(R)) {
      R <- R + w
      K <- K - 1
    }
  }
  repeat {
    x1 <- L + runif(1) * (R - L)
    if (y < f(x1)) 
      return(list(x = x1, x0 = x, l = L, r = R, y = y, nEvaluations = nEvaluations))
    if (x1 < x) 
      L <- x1
    else R <- x1
  }
}

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
samples <- list(beta = matrix(0.0, nrow = n_samples, ncol = length(beta)),
                psi = numeric(n_samples), g = numeric(n_samples),
                time = numeric(1), w = numeric(1), nEval = numeric(n_samples),
                x0 = numeric(n_samples), y = numeric(n_samples),
                l = numeric(n_samples), r = numeric(n_samples))
# setting up chain storage
chainSamples <- vector('list', length = nChains)
chainSamples <- lapply(chainSamples, \(list) samples)

xx <- seq(from = 0.01, to = 100, length.out = 1000)

f <- \(g) {
  beta_cov <- g / psi * inv_XtX
  exp(dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf))
}
yy <- sapply(xx, f)
plot(xx,yy)

w <- c(30, 35, 40, 45, 50)

unlink('explorationPlots/*', recursive = TRUE)

sapply(w, FUN = \(w) dir.create(paste0('explorationPlots/w',w)))

# collecting the samples
output <- foreach( chain = seq_along(chainSamples) ) %do% {
# for( chain in seq_along(chainSamples) ) {
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
      # updating g
      temp <- slice_sampler_stepping_out(g, \(g) {
        beta_cov <- g / psi * inv_XtX
        dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf)
      }, w = w[chain], log = TRUE, max = Inf)
      # plotting the update process
      if(i == 1) { par(mfrow = c(4,5)); pdf(file = paste0('explorationPlots/w',w[chain],'/',i,'iter.pdf')) }
      if(i %% 1000 == 0 | i %in% c(1,2,3,4,5)) {
        f <- \(g) {
          beta_cov <- g / psi * inv_XtX
          exp(dmvnorm(beta, beta_0, beta_cov, log = TRUE) + ifelse(0 < g && g < g_max, 0, -Inf))
        }
        yy <- sapply(xx, f)
        expy <- exp(temp$y)
        plot(xx,yy, main = paste0('g=',round(g),', iter=',i),cex = 0.5)
        points(x = temp$x, y = expy, col = 'red') # new point
        points(x = temp$x0, y = expy, col = 'dodgerblue') # original point
        points(x = c(temp$l,temp$r), y = c(expy,expy), col = 'black', pch = 3)
        # readline(prompt = 'continue?')
      }
      if(i == n_samples) dev.off()
      # saving all the variables in the plotting process
      g <- temp$x
      chainSamples[[chain]]$beta[i,] <- beta
      chainSamples[[chain]]$psi[i] <- psi
      chainSamples[[chain]]$g[i] <- g
      chainSamples[[chain]]$nEval[i] <- temp$nEvaluations
      chainSamples[[chain]]$x0[i] <- temp$x0
      chainSamples[[chain]]$y[i] <- temp$y
      chainSamples[[chain]]$l[i] <- temp$l
      chainSamples[[chain]]$r[i] <- temp$r
    }
  })
  chainSamples[[chain]]$time <- time['user.self']
  chainSamples[[chain]]$w <- w[chain]
}


# evaluation of samples

# tst <- sapply(1:ncol(gSamples) , FUN = \(i) {
#   col <- gSamples[,i]
#   effSize <- round(coda::effectiveSize(col))
#   sampPSec <- round(effSize/time[i])
#   plot(density(col), main = paste0('w=',w[i],', SampPSec=',sampPSec))
#   fiveSummary <- round(summary(col))
#   plot(col, type = 'l', main = paste0('w=',w[i],', min=',fiveSummary[1],', mean=',fiveSummary[4],', max=',fiveSummary[5]))
#   data.frame(w = w[i], sampPSec = sampPSec, est = round(mean(col),2), lwrCred = round(quantile(col, probs = 0.025)), uprCred = round(quantile(col, probs = 0.975)) )
# })
# 
# t(tst)


saveRDS(chainSamples, file = 'data/steppingOutExploration.rds')

print('finished')
