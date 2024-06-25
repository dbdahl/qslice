
lapx_g <- function(p, psi, Q, a, sc_adj = 1.0, logG = FALSE) {

  if (isTRUE(logG)) {

    A <- a + p - 2
    psiQ <- psi * Q
    negB <- (psiQ - p + 2)

    quad_sol <- (negB + c(-1, 1) * sqrt( negB^2 + 4*A*psiQ )) / (2 * A)
    eyhat <- max(quad_sol)
    sc <- sc_adj / sqrt( 0.5 * ( a * eyhat / (1.0 + eyhat)^2 + psiQ / eyhat ) )

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
  pseudo_list(family = "t",
              params = list(loc = tmp$loc, sc = tmp$sc,
                            degf = degf),
              lb = lb, ub = ub, name = 'Laplace')
}
