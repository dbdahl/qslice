prior <- list()

## beta ~ mvn(beta_0, g / psi * inv(XtX)), using the covariance parametrization
prior$beta_0 <- rep(0.0, dat$p)

## psi ~ gamma(a0, b0) => sig2 ~ ig(a_0, b_0)
prior$n_0 <- 5
prior$s_0 <- 0.4

prior$a_0 <- prior$n_0 / 2.0
prior$b_0 <- prior$n_0 * prior$s_0^2 / 2.0

## p(g) propto (1 + g)^(-a/2) I(0, g_max) ; g ~ uniform(0, g_max) if a = 0
prior$g_max <- 3 * dat$p^2
prior$a_g <- 3

prior
