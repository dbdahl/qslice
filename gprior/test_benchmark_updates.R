
update_beta <- function(g, psi, 
                        inv_XtX, 
                        beta_mle, 
                        beta_0) {
  require("mvtnorm")
  q <- g / (1.0 + g)
  beta_cov <- q * inv_XtX / psi
  beta_mean <- q * beta_mle + (1.0 - q) * beta_0
  mvtnorm::rmvnorm(1, beta_mean, beta_cov) |> drop()
}

# need Li Li' = inv_XtX
# XtX = inv(Li Li') = inv(Li') inv(Li) = inv(Li)' inv(Li)
# how do we get Li from R (upper triangular) where R'R = XtX?
# inv(XtX) = inv(R'R) = inv(R) inv(R)' ==> inv(R) = Li'? NO
# can we get it from LL' = XtX?
# inv(XtX) = inv(LL') = inv(L)' inv(L) ==> inv(R) = inv(L)'? YES, but these are upper triangular (doesn't matter)


update_beta2 <- function(g, psi, 
                         inv_chol_XtX, 
                         beta_mle, 
                         beta_0) {
  q <- g / (1.0 + g)
  beta_mean <- q * beta_mle + (1.0 - q) * beta_0
  z <- rnorm(length(beta_0))
  out <- sqrt(q/psi) * inv_chol_XtX %*% z + beta_mean
  drop(out)
}

update_beta3 <- function(g, psi, 
                         chol_XtX, 
                         beta_mle, 
                         beta_0) {
  q <- g / (1.0 + g)
  beta_mean <- q * beta_mle + (1.0 - q) * beta_0
  z <- rnorm(length(beta_0))
  out <- sqrt(q/psi) * backsolve(chol_XtX, z) + beta_mean
  drop(out)
}


g <- 2
psi <- 0.3
p <- 15
n <- 20
X <- matrix(rnorm(n*p), nrow = n)
XtX <- crossprod(X)
chol_XtX <- chol(XtX)
inv_XtX <- chol2inv(chol_XtX)
inv_chol_XtX <- backsolve(chol_XtX, diag(p))
beta_0 <- rep(1, p)
beta_mle <- rnorm(p, 1, 0.1)

library("microbenchmark")

n_test <- 100e3
Y1 <- replicate(n_test, update_beta(g, psi, inv_XtX, beta_mle, beta_0)) |> t()
Y2 <- replicate(n_test, update_beta2(g, psi, inv_chol_XtX, beta_mle, beta_0)) |> t()
Y3 <- replicate(n_test, update_beta3(g, psi, chol_XtX, beta_mle, beta_0)) |> t()

colMeans(Y1)
colMeans(Y1) - colMeans(Y2); max(abs(colMeans(Y1) - colMeans(Y2)))
colMeans(Y1) - colMeans(Y3); max(abs(colMeans(Y1) - colMeans(Y3)))

cov(Y1)
cov(Y1) - cov(Y2); max(abs(cov(Y1) - cov(Y2)))
cov(Y1) - cov(Y3); max(abs(cov(Y1) - cov(Y3)))

microbenchmark(
  update_beta(g, psi, inv_XtX, beta_mle, beta_0),
  update_beta2(g, psi, inv_chol_XtX, beta_mle, beta_0),
  update_beta3(g, psi, chol_XtX, beta_mle, beta_0)
)

microbenchmark(
  update_beta2(g, psi, inv_chol_XtX, beta_mle, beta_0),
  update_beta3(g, psi, chol_XtX, beta_mle, beta_0)
)


cp1 <- function(X, bd) {
  crossprod(X %*% bd) |> drop()
}

cp2 <- function(XtX, bd) {
  bd %*% XtX %*% bd |> drop()
}

bd <- 1:p

microbenchmark(
  cp1(X, bd),
  cp2(XtX, bd),
  times = 10e3
)
