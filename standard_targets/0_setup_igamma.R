## Set up for the Inverse gamma
## author: Sam Johnnson

source("functions/invgamma.R")

gshape <- 2.0

truth = list(d = function(x) {dinvgamma(x, shape = gshape, scale = 1.0)},
             ld = function(x) {dinvgamma(x, shape = gshape, scale = 1.0, log = TRUE)},
             p = function(x) {pinvgamma(x, shape = gshape, scale = 1.0)},
             q = function(u) {qinvgamma(u, shape = gshape, scale = 1.0)},
             t = paste0("inv-gamma(shape = ", gshape, ")"),
             lb = 0.0,
             ub = Inf)
