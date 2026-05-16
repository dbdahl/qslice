
# data_use <- "mtcars"
# data_use <- "diabetes"
data_use <- "riboflavin"

source("0_data.R")
source("0_prior.R")


## STAN
library("rstan")

data_stan <- list(y = c(dat$y), X = as.matrix(dat$X), n = dat$n, p = dat$p,
                  n0 = prior$n0, s02 = prior$s02)

params <- c("beta", "llam2", "tau2", "ltau2", "sig2")

# n_cores <- parallel::detectCores()   # Count how many cores are available
n_cores <- 3
options(mc.cores = n_cores)

fit_stan <- stan(model_code = readLines("hs.stan"),
                 data = data_stan, pars = params,
                 iter = 5000, warmup = 500, thin = 5, chains = n_cores)

sim_stan <- extract(fit_stan)

# summary(fit_stan)
plot(fit_stan)

plot(fit_stan, pars = "llam2")

rstan::traceplot(fit_stan, pars = "ltau2")
rstan::traceplot(fit_stan, pars = "beta[73]")


(nsim <- length(sim_stan$tau))

summary(fit_stan)$s
