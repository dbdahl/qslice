
lapx_g <- function(p, psi, Q, a, sc_adj = 1.0, logG = FALSE) {

  if (isTRUE(logG)) {

    A <- a + p - 2
    psiQ <- psi * Q
    negB <- (psiQ - p + 2)

    # if (a > 0) {
    quad_sol <- (negB + c(-1, 1) * sqrt( negB^2 + 4*A*psiQ )) / (2 * A)
    eyhat <- max(quad_sol)
    sc <- sc_adj / sqrt( 0.5 * ( a * eyhat / (1.0 + eyhat)^2 + psiQ / eyhat ) )
    # } else {
    # xhat <- psi * Q / p # when a = 0
    # sc <- sc_adj * xhat * sqrt(2.0 / p) # when a = 0
    # }

    out <- list(loc = log(eyhat), sc = sc)

  } else {

    A <- a + p
    psiQ <- psi * Q
    negB <- (psiQ - p)

    if (a > 0) {
      quad_sol <- (negB + c(-1, 1) * sqrt( negB^2 + 4*A*psiQ )) / (2 * A)
      xhat <- max(quad_sol)
      sc <- sc_adj / sqrt( psiQ / xhat^3 - p / (2 * xhat^2) - a / (2 * (1 + xhat)^2) )
    } else {
      xhat <- psi * Q / p # when a = 0
      sc <- sc_adj * xhat * sqrt(2.0 / p) # when a = 0
    }

    out <- list(loc = xhat, sc = sc)

  }

  out
}

lapxt_g <- function(p, psi, Q, a, sc_adj = 1.0, degf = 1.0,
                    lb = -Inf, ub = Inf, logG = FALSE) {
  tmp <- lapx_g(p = p, psi = psi, Q = Q, a = a, sc_adj = sc_adj, logG = logG)
  pseudo_t_list(loc = tmp$loc, sc = tmp$sc,
                degf = degf, lb = lb, ub = ub, name = 'Laplace')
}


### tests
lf <- function(x, p, psi, Q, a, logG = FALSE) {
  out <- 0.0
  if (isTRUE(logG)) {
    out <- out + x # right now x = log(g); this is pre-adding the Jacobian
    x <- exp(x) # then x becomes x
  }
  out <- out - p*log(x)/2 - a*log1p(x)/2 - psi*Q/(2*x)
}

# p <- 10
# psi <- 8
# Q <- 20
# a <- 3
# sc_adj <- 1.2
# logG <- FALSE
# degf <- 5
#
# xlim <- c(1, 200)
# # xlim <- c(1, 1e3)
# # xlim <- c(1e10, 1e12)
# if (logG) { xlim <- log(xlim)}
#
# (lapx <- lapx_g(p, psi, Q, a, sc_adj = sc_adj, logG = logG))
# (lapx0 <- lapx_g(p + a, psi, Q, 0, sc_adj = sc_adj, logG = logG))
#
# # (lapx <- lapproxt(lf = lf, init = 10.0, sc_adj = 1, lb = 0.0, p = p, psi = psi, Q = Q))
# # (lapx <- lapproxt(lf = lf, init = 10.0, sc_adj = 2, lb = 0.0, p = p, psi = psi, Q = Q))
#
# curve(exp(lf(x, p, psi, Q, a, logG) - lf(lapx$loc, p, psi, Q, a, logG)), from = xlim[1], to = xlim[2], n = 1e3,
#       main = paste0("p = ", p, "  psi = ", psi, "  Q = ", Q))
# curve(exp(lf(x, p + a, psi, Q, 0, logG) - lf(lapx0$loc, p + a, psi, Q, 0, logG)), from = xlim[1], to = xlim[2], n = 1e3, add = TRUE, lty = 2)
# abline(v = lapx$loc, lty=2)
# curve(dt((x - lapx$loc)/lapx$sc, df = degf) / dt(0, df = degf), from = xlim[1], to = xlim[2], n = 1e3,
#       add = T, col = "red")
#
# lapxt_g(p, psi, Q, a, sc_adj = 1.0, lb = 0.0, logG = logG)
