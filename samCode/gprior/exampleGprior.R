library(mvtnorm)
library(cucumber)

data(mtcars)

attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)

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

n_samples <- 10000
samples <- list(beta = matrix(0.0, nrow = n_samples, ncol = length(beta)), psi = numeric(n_samples), g = numeric(n_samples))
for (i in seq_len(n_samples)) {
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
  }, w = 100, log = TRUE)$x
  samples$beta[i,] <- beta
  samples$psi[i] <- psi
  samples$g[i] <- g
}

plot(density(samples$g))
acf(samples$g)
coda::effectiveSize(samples$g)


