gshape <- 2.5

truth <- list(d = function(x) {dgamma(exp(x), gshape) * exp(x)},
              ld = function(x) {dgamma(exp(x), gshape, log = TRUE) + x},
              # dld = function(x) {-x},
              p = function(x) {pgamma(exp(x), shape = gshape)},
              q = function(u) {log(qgamma(u, shape = gshape))},
              lb = -Inf, ub = Inf,
              t = paste("gamma_log(shape = ", gshape, ")")
              )
