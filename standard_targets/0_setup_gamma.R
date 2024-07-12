## Set up for the cauch transform
## author: Sam Johnnson

gshape <- 2.5

truth = list(d = function(x) {dgamma(x, shape = gshape, scale = 1.0) * (x > 0.0)},
             ld = function(x) {dgamma(x, shape = gshape, scale = 1.0,  log=TRUE) + ifelse(x > 0.0, 0.0, -Inf)},
             # dld = function(x) {(gshape - 1.0)/x - 1.0},
             p = function(x) {pgamma(x, shape = gshape, scale = 1.0)},
             q = function(u) {qgamma(u, shape = gshape, scale = 1.0)},
             t = paste0("gamma(shape = ", gshape, ")"),
             lb = 0.0,
             ub = Inf)

xlim_range <- c(-4, 15)
ylim_range <- c(0, 0.42)
