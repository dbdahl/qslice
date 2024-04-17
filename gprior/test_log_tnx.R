
f <- function(x, p, psi, Q, a) {
  x^(-0.5*p) * (1 + x)^(-0.5*a) * exp(-0.5*psi*Q/x)
}

p <- 10
psi <- 10
Q <- 25
a <- 3

curve(f(x, p, psi, Q, a), from = 0, to = 200, n = 1e3)

flx <- function(x, p, psi, Q, a) f(exp(x), p, psi, Q, a) * exp(x)

curve(flx(x, p, psi, Q, a), from = 0, to = 7, n = 1e3)
