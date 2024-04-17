
if (type == "AR") {

  (Cor <- outer(1:K, 1:K, FUN = function(i, j) rho^abs(i-j)))
  sig <- 1
  Sig <- sig * Cor
  SigL <- chol(Sig) |> t()

  Y <- SigL %*% matrix(rnorm(n*K), nrow = K) |> t()
  ltarget <- function(x) -0.5 * drop(crossprod(forwardsolve(SigL, x)))

} else if (type == "indep") {

  sig2 <- rgamma(K, shape = 2, scale = 5)
  Sig <- diag(sig2)
  SigL <- chol(Sig)

  Y <- SigL %*% matrix(rnorm(n*K), nrow = K) |> t()
  ltarget <- function(x) -0.5 * sum(x^2 / sig2)

} else if (type == "corrBeta") {
  stopifnot(K == 2)

  Y <- matrix(NA, nrow = n, ncol = K)
  Y[,1] <- rbeta(n, 2, 3)
  Y[,2] <- rbeta(n, abs((Y[,1] + 0.5*rho^2)^(6*rho)/(0.15*rho^2) + 2), 4)

  ltarget <- function(x) {
    if (any(x < 0.0)) {
      out <- -Inf
    } else {
      out <- dbeta(x[1], 2, 3, log = TRUE) +
             dbeta(x[2], abs((x[1] + 0.5*rho^2)^(6*rho)/(0.15*rho^2) + 2), 4, log = TRUE)
    }
    out
  }

  pseudo_control <- list(

    loc_fn = function(x) {
      if (is.null(x)) {
        out <- 0.5
      } else {
        out <- 0.13 * log(x[length(x)]) + 1.0
      }
      out
    },

    sc_fn = function(x) {
      if (is.null(x)) {
        out <- 0.3
      } else {
        out <- 0.6 / (10.0*x[length(x)] + 1.0)
      }
      out
    },

    pseu_init = pseudo_t_list(loc = 0.5, sc = 0.3, degf = 5.0,
                              lb = 0.0, ub = 1.0),

    degf = 5.0,
    lb = rep(0.0, K),
    ub = rep(1.0, K)
  )

} else if (type == "funnel") {

  Y <- matrix(NA, nrow = n, ncol = K)
  Y[,1] <- rnorm(n, 0.0, 3.0)

  for (i in 1:n) {
    Y[i, 2:K] <- rnorm(K-1, 0.0, exp(0.5*Y[i,1]))
  }

  ltarget <- function(x) {
    dnorm(x[1], 0.0, 3.0, log = TRUE) +
      sum(dnorm(x[2:K], 0.0, exp(0.5*x[1]), log = TRUE))
  }

  pseudo_control <- list(

    loc_fn = function(x) {
      0.0
    },

    sc_fn = function(x) {
      if (is.null(x)) {
        out <- 3.0
      } else {
        out <- exp(0.5*x[1])
      }
      out
    },

    pseu_init = pseudo_t_list(loc = 0.0, sc = 3.0, degf = 5.0,
                              lb = -Inf, ub = Inf),

    degf = 5.0,
    lb = rep(-Inf, K),
    ub = rep(Inf, K)
  )

} else if (type == "funnel_half") {

  Y <- matrix(NA, nrow = n, ncol = K)
  Y[,1] <- rnorm(n, 0.0, 3.0)

  for (i in 1:n) {
    Y[i, 2:K] <- rnorm(K-1, 0.0, exp(0.5*Y[i,1])) |> abs()
  }

  ltarget <- function(x) {
    if(any(x[-1] <= 0.0)) {
      out <- -Inf
    } else {
      out <- dnorm(x[1], 0.0, 3.0, log = TRUE) +
        sum(dnorm(x[2:K], 0.0, exp(0.5*x[1]), log = TRUE))
    }
    out
  }

  pseudo_control <- list(

    loc_fn = function(x) {
      0.0
    },

    sc_fn = function(x) {
      if (is.null(x)) {
        out <- 3.0
      } else {
        out <- exp(0.5*x[1])
      }
      out
    },

    pseu_init = pseudo_t_list(loc = 0.0, sc = 3.0, degf = 5.0,
                              lb = -Inf, ub = Inf),

    degf = 5.0,
    lb = c(-Inf, rep(0.0, K-1)),
    ub = rep(Inf, K)
  )
}

str(Y)
cov(Y)
plot(Y[,1], Y[,2])

(mu_hat <- colMeans(Y))
(Sig_hat <- cov(Y))
(SigL_hat <- t(chol(Sig_hat)))
