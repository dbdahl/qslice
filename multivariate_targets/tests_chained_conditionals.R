
## these go with the teardrop beta with rho = 0.25
loc_fn <- function(x) {
  if (is.null(x)) {
    out <- 0.5
  } else {
    out <- 0.13 * log(x[length(x)]) + 1.0
  }
  out
}

sc_fn <- function(x) {
  if (is.null(x)) {
    out <- 0.3
  } else {
    out <- 0.6 / (10.0*x[length(x)] + 1.0)
  }
  out
}


ps_init <- pseudo_t_list(loc = loc_fn(NULL), sc = sc_fn(NULL), degf = degf, 
                         lb = 0.0, ub = 1.0)

# x <- c(0.1, 0.6)
# x <- c(0.4, 0.6)
# x <- c(0.7, 0.6)
# x <- c(0.9, 0.6)
# 
# ps_seq <- pseudo_t_condseq(x, init_t = ps_init, 
#                  loc_fn = loc_fn, sc_fn = sc_fn, degf = 5, 
#                  lb = rep(0.0, length(x)), ub = rep(1.0, length(x)))

n_sim <- 10e3
Ystar <- matrix(NA, nrow = n_sim, ncol = K)
U <- matrix(runif(n_sim*K), nrow = n_sim)

for (i in 1:n_sim) { # get x from u (i.e., simulate from pseudo)
  
  # Ystar[i, 1] <- ps_init$q(U[i,1])
  # ps_seq <- pseudo_t_condseq(Ystar[i, 1], init_t = ps_init, 
  #                            loc_fn = loc_fn, sc_fn = sc_fn, degf = 5, 
  #                            lb = rep(0.0, K), ub = rep(1.0, K))
  # Ystar[i, 2] <- ps_seq[[2]]$q(U[i,2])
  
  Ystar[i,] <- pseudo_t_condseq_XfromU(U[i,], init_t = ps_init, 
                                       loc_fn = loc_fn, sc_fn = sc_fn, 
                                       degf = 5, 
                                       lb = rep(0.0, K), ub = rep(1.0, K))$x
  cat(i, "\r")
}

plot(Y[,1], Y[,2]) # generated from 0_targets.R
points(Ystar[,1], Ystar[,2], col = "red", cex = 0.6)
summary(Ystar)

Ustar <- matrix(NA, nrow = n_sim, ncol = K)
for (i in 1:n_sim) { # get u from x
  ps_seq <- pseudo_t_condseq(Ystar[i, ], init_t = ps_init, 
                             loc_fn = loc_fn, sc_fn = sc_fn, degf = 5, 
                             lb = rep(0.0, K), ub = rep(1.0, K))
  Ustar[i, ] <- sapply(1:K, function(k) ps_seq[[k]]$p(Ystar[i,k]))
  cat(i, "\r")
}

plot(Ustar)
max(abs(Ustar - U)) # recovered!

