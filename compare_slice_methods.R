rm(list=ls())

library(ninja)
source("functions.R")

# Log of density of the target distribution (select one)

f <- function(x) {
    counter <<- counter + 1
    log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2))
}

f <- function(x) {
  counter <<- counter + 1
  log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=20,sd=1))
}

f <- function(x) {
  counter <<- counter + 1
  dgamma(x, shape=2.5, log=TRUE)
}

f <- function(x) {
  counter <<- counter + 1
  ifelse( x <= 0, -Inf, dt(x, df=1.0, log=TRUE) + log(2.0))
}

f <- function(x) {
  counter <<- counter + 1
  dt(x, df=1.0, log=TRUE)
}

f <- function(x) {
  counter <<- counter + 1
  dbeta(x, shape1=0.2, shape2=0.8, log=TRUE)
}


## Pseudo prior (select one)
logPseudoPDF <- function(x) dnorm(x, mean=4.0, sd=15.0, log=TRUE)
pseudoInvCDF <- function(u) qnorm(u, mean=4.0, sd=15.0)

logPseudoPDF <- function(x) dnorm(x, mean=4.8, sd=3.0, log=TRUE)
pseudoInvCDF <- function(u) qnorm(u, mean=4.8, sd=3.0)

logPseudoPDF <- function(x) dnorm(x, mean=0, sd=1.5, log=TRUE)
pseudoInvCDF <- function(u) qnorm(u, mean=0, sd=1.5)

logPseudoPDF <- function(x) dt((x-4.8)/3, df=5, log=TRUE) - log(3)
pseudoInvCDF <- function(u) 4.8 + 3*qt(u, df=5)

logPseudoPDF <- function(x) dgamma(x, shape=2.5, log=TRUE)
pseudoInvCDF <- function(u) qgamma(u, shape=2.5)

logPseudoPDF <- function(x) dexp(x, log=TRUE)
pseudoInvCDF <- function(u) qexp(u)

logPseudoPDF <- function(x) if ( x <=0 ) -Inf else dt(x, df=1.0, log=TRUE) + log(2.0)
pseudoInvCDF <- function(u) qt((u + 1.0)/2.0, df=1) # half Cauchy

logPseudoPDF <- function(x) dbeta(x, shape1=0.5, shape2=0.5, log=TRUE)
pseudoInvCDF <- function(u) qbeta(u, shape1=0.5, shape2=0.5)




## Samples with stepping-out procedure
counter <- 0
draws0 <- numeric(50000)
draws0[1] <- 0.5
time0 <- system.time({
    for ( i in seq_along(draws0)[-1] ) {
      draws0[i] <- slice_sampler(draws0[i-1], f, w=20, max=Inf, log=TRUE)
    }
})
counter0 <- counter


## Samples with proposed method
counter <- 0
f_compos <- function(u) {
  x <- pseudoInvCDF(u)
  f(x) / pseudoPDF(x)
}

draws1 <- numeric(50000)
time1 <- system.time({
  u <- 0.5
  fx <- f_compos(u)
  if ( fx <= 0 ) stop("Oops, the starting value is not in support.")
  current <- list(x=u, fx=fx)
  for ( i in seq_along(draws1) ) {
    current <- slice_sampler_shrinkage(current, f_compos)
    draws1[i] <- pseudoInvCDF(current$x)
  }
})
counter1 <- counter


## Samples with latent slice method (Li & Walker, 2020)
counter <- 0
draws2 <- numeric(50000)
time2 <- system.time({
  x <- 0.1
  s <- 1.0
  lambda <- 0.1
  fx <- f(x)
  if ( fx <= 0 ) stop("Oops, the starting value is not in support.")
  current <- list(x=x, s=s, fx=fx)
  for ( i in seq_along(draws0) ) {
    current <- slice_latent(current, f, lambda)
    draws2[i] <- current$x
  }
})
counter2 <- counter




## check
# Stepping out
hist(draws0[-c(1:10000)], freq=FALSE, breaks=50)
lines(density(draws0), col="red")
curve(exp(f(x)), -3, 25, add=TRUE, col="black", n=1000)
(B0 <- coda::effectiveSize(draws0))   # Effective number of independent samples
B0 / time0['user.self']             # Effective number of independent samples per CPU second
counter0                           # Number of function evaluations

x_test = 0.7

integrate(function(x) exp(f(x)), -Inf, x_test)
mean(draws0 < x_test)
(err0 <- abs( integrate(function(x) exp(f(x)), -Inf, x_test)$value - mean(draws0 < x_test) ))

(acf0 <- acf(draws0)$acf[2,1,1])


# Transform
hist(draws1[-c(1:10000)], freq=FALSE, breaks=50)
lines(density(draws1), col="red")
curve(f, -3, 25, add=TRUE, col="black", n=1000)
curve(pseudoPDF, -3, 25, add=TRUE, col="blue", lty=2)
(B1 <- coda::effectiveSize(draws1))  # Effective number of independent samples
B1 / time1['user.self']             # Effective number of independent samples per CPU second
counter1                            # Number of function evaluations

integrate(f, -Inf, x_test)
mean(draws1 < x_test)
(err1 <- abs( integrate(f, -Inf, x_test)$value - mean(draws1 < x_test) ))

(acf1 <- acf(draws1)$acf[2,1,1])


# Latent
hist(draws2[-c(1:10000)], freq=FALSE, breaks=50)
lines(density(draws2), col="red")
curve(f, -3, 25, add=TRUE, col="black", n=1000)
(B2 <- coda::effectiveSize(draws2))  # Effective number of independent samples
B2 / time2['user.self']             # Effective number of independent samples per CPU second
counter2                            # Number of function evaluations

integrate(f, -Inf, x_test)
mean(draws2 < x_test)
(err2 <- abs( integrate(f, -Inf, x_test)$value - mean(draws2 < x_test) ))

(acf2 <- acf(draws2)$acf[2,1,1])


## compare
err0 # stepping out
err1 # transform
err2 # latent

acf0
acf1
acf2

(B0 <- coda::effectiveSize(draws0))   # Effective number of independent samples
(B1 <- coda::effectiveSize(draws1))   # Effective number of independent samples
(B2 <- coda::effectiveSize(draws2))   # Effective number of independent samples

B0 / time0['user.self']
B1 / time1['user.self']
B2 / time2['user.self']

counter0 # stepping out
counter1 # transform
counter2 # latent
counter0 / counter1

