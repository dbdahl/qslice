rm(list = ls()); dev.off()
library("cucumber")
library("coda")
source("0_samplers.R")

type <- "AR"; rho <- 0.8
type <- "indep"
type <- "corrBeta"; rho <- 0.25 # 1.0 or 0.25
type <- "funnel"
type <- "funnel_half"

K <- 2
n <- 15e3
source("0_targets.R")
ls()


n_iter <- 10e3

degf <- 1



sc_adj <- 1.3
(SigL_use <- sc_adj * SigL_hat)
(SigL_use <- diag(sc_adj * sqrt(diag(Sig_hat)))) # the target uses SigL, which is also defined in the sampler; correctness depends on lexical scoping


time <- system.time({
  mcmc_out <- gess_sampler(ltarget, n_iter = n_iter, x0 = runif(K),
                           loc = mu_hat, ScL = SigL_use, degf = degf, is_chol = TRUE)
})

(es <- effectiveSize(as.mcmc(mcmc_out$draws)))
mean(es)
mean(es) / time["user.self"]

mcmc_out$n_eval / n_iter

plot(mcmc_out$draws[,1], mcmc_out$draws[,2])
points(Y[,1], Y[,2], col = "blue", cex = 0.5)

cov(mcmc_out$draws)
colMeans(mcmc_out$draws)


## chain of indep pseudo
ps <- lapply(1:K, function(k) {
  list(ld = function(x) dt((x - mu_hat[k])/SigL_use[k,k], df = degf, log = TRUE),
       p = function(x) pt((x - mu_hat[k])/SigL_use[k,k], df = degf),
       q = function(u) mu_hat[k] + SigL_use[k,k] * qt(u, df = degf) )
})

time <- system.time({
  mcmc_out <- Qslice_sampler(ltarget, n_iter = n_iter, x0 = runif(K), pseudo = ps)
})

time <- system.time({
  mcmc_out <- Qslice_seq_sampler(ltarget, n_iter = n_iter, x0 = runif(K),
                                 pseudo_control = pseudo_control)
})

(es <- effectiveSize(as.mcmc(mcmc_out$draws)))
mean(es)
mean(es) / time["user.self"]

mcmc_out$n_eval / n_iter
sapply(1:K, function (k) auc(u = mcmc_out$draws_u[,k]))

hist(mcmc_out$draws_u[,2])

plot(mcmc_out$draws[,1], mcmc_out$draws[,2])
points(Y[,1], Y[,2], col = "blue", cex = 0.5)
cov(mcmc_out$draws)
colMeans(mcmc_out$draws)
