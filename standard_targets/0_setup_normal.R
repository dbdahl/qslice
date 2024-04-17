## authors: Sam Johnnson, Matt Heiner

truth <- list(d = function(x) dnorm(x),
             ld = function(x) dnorm(x, log=TRUE),
             # dld = function(x) {-x},
             p = function(x) pnorm(x),
             q = function(u) qnorm(u),
             lb = -Inf, ub = Inf,
             t = "normal(0,1)"
             )
truth$dld <- function(x) numDeriv::grad(truth$ld, x=x)

xlim_range <- c(-4, 4)
ylim_range <- c(0, 0.42)
