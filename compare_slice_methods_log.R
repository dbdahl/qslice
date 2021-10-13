rm(list=ls())

library(cucumber)
source("functions_log.R")

# Log of density of the target distribution (select one)

lf <- function(x) {
    log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2))
}

lf <- function(x) {
    log(0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=20,sd=1))
}

lf <- function(x) {
    dgamma(x, shape=2.5, log=TRUE)
}

lf <- function(x) {
    ifelse( x < 0, -Inf, dt(x, df=1.0, log=TRUE) + log(2.0))
}

lf <- function(x) {
    dt(x, df=1.0, log=TRUE)
}

lf <- function(x) {
    dbeta(x, shape1=0.2, shape2=0.8, log=TRUE)
}

f <- function(x) exp(lf(x))

## Samples with stepping-out procedure (Neal 2003)
counter <- 0
draws <- numeric(50000)
draws[1] <- 0.5
time <- system.time({
    for ( i in seq.int(2,length(draws)) ) {
        out <- slice_sampler_stepping_out(draws[i-1], lf, w=10, max=Inf, log=TRUE)
        draws[i] <- out$x
        counter <- counter + out$nEvaluations
    }
})
counter
plot(density(draws))
coda::effectiveSize(draws) / time['user.self']

## Samples with latent procedure (Li and Walker 2020)
counter <- 0
draws <- latents <- numeric(50000)
draws[1] <- 0.5
latents[1] <- 0.3
time <- system.time({
    for ( i in seq.int(2,length(draws)) ) {
        out <- slice_sampler_latent(draws[i-1], latents[i-1], lf, rate=0.03)
        draws[i] <- out$x
        latents[i] <- out$s
        counter <- counter + out$nEvaluations
    }
})
counter
plot(density(draws))
coda::effectiveSize(draws) / time['user.self']



######
###### The code below is not currently tested / functional
######


# Setup for elliptical slice sampler
ellip_ctr <- 0.0
ellip_sc <- 1.0
ellip_degf <- 10.0

ellip_invTform <- function(z) z # identity transformation
ellip_logJacobian <- function(z) 0.0

ellip_invTform <- function(z) exp(z) # log transformation
ellip_logJacobian <- function(z) z

ellip_invTform <- function(z) 1.0 / (1.0 + exp(-z)) # logit transformation
ellip_logJacobian <- function(z) {
    enz <- exp(-z)
    -z - 2.0*log(1.0 + enz)
}

ellip_logtarget <- function(z) lf( ellip_invTform(z) ) + ellip_logJacobian(z)


## Pseudo prior (select one)
pseudoLogPDF <- function(x) dnorm(x, mean=4.0, sd=15.0, log=TRUE)
pseudoInvCDF <- function(u) qnorm(u, mean=4.0, sd=15.0)

pseudoLogPDF <- function(x) dt((x-ellip_ctr)/ellip_sc, df=ellip_degf, log=TRUE) - log(ellip_sc)
pseudoInvCDF <- function(u) ellip_sc * qt(u, df=ellip_degf) + ellip_ctr

pseudoLogPDF <- function(x) dgamma(x, shape=2.5, log=TRUE)
pseudoInvCDF <- function(u) qgamma(u, shape=2.5)

pseudoLogPDF <- function(x) dexp(x, log=TRUE)
pseudoInvCDF <- function(u) qexp(u)

pseudoLogPDF <- function(x) if ( x < 0 ) -Inf else dt(x, df=1.0, log=TRUE) + log(2.0)
pseudoInvCDF <- function(u) qt((u + 1.0)/2.0, df=1) # half Cauchy

pseudoLogPDF <- function(x) dbeta(x, shape1=0.5, shape2=0.5, log=TRUE)
pseudoInvCDF <- function(u) qbeta(u, shape1=0.5, shape2=0.5)




## Samples with proposed method
counter <- 0
lf_compos <- function(u) {
    x <- pseudoInvCDF(u)
    lf(x) - pseudoLogPDF(x)
}
draws1 <- numeric(50000)
time1 <- system.time({
    u <- 0.5
    lfx <- lf_compos(u)
    if ( lfx == -Inf ) stop("Oops, the starting value is not in support.")
    current <- list(x=u, lfx=lfx)
    for ( i in seq_along(draws1) ) {
        current <- slice_sampler_shrinkage(current, lf_compos)
        draws1[i] <- pseudoInvCDF(current$x)
    }
})
counter1 <- counter

counter <- 0
draws10 <- numeric(50000)
draws10[1] <- 0.5
time10 <- system.time({
    for ( i in seq_along(draws1)[-1] ) {
        draws10[i] <- slice_sampler_transform(draws10[i-1], lf, pseudoLogPDF, pseudoInvCDF)
    }
})
counter10 <- counter


## Samples with latent slice method (Li & Walker, 2020)
counter <- 0
draws2 <- numeric(50000)
time2 <- system.time({
    x <- 0.1
    s <- 1.0
    lambda <- 0.1
    lfx <- lf(x)
    if ( lfx == -Inf ) stop("Oops, the starting value is not in support.")
    current <- list(x=x, s=s, lfx=lfx)
    for ( i in seq_along(draws1) ) {
        current <- slice_latent(current, lf, lambda)
        draws2[i] <- current$x
    }
})
counter2 <- counter


## Samples with elliptical slice sampler
draws3 <- numeric(50000)
time3 <- system.time({
    z <- 0.1
    lfz <- ellip_logtarget(z)
    if ( lfz == -Inf ) stop("Oops, the starting value is not in support.")
    current <- list(x=z, lfx=lfz)
    for ( i in seq_along(draws3) ) {
        current <- slice_ellipse_generalized(current, ellip_logtarget, ellip_ctr, ellip_sc, ellip_degf)
        draws3[i] <- ellip_invTform(current$x)
    }
})
counter3 <- counter

## version of elliptical sampler with no transform (considerably faster than identity transform)
# draws2 <- numeric(50000)
# time2 <- system.time({
#   z <- 0.1
#   lfz <- lf(z)
#   if ( lfz == -Inf ) stop("Oops, the starting value is not in support.")
#   current <- list(x=z, lfx=lfz)
#   for ( i in seq_along(draws2) ) {
#     current <- slice_ellipse_generalized(current, lf, ellip_ctr, ellip_sc, ellip_degf)
#     draws2[i] <- current$x
#   }
# })
# counter2 <- counter
# counter <- 0




## check
# Stepping out
hist(draws0[-c(1:10000)], freq=FALSE, breaks=50)
lines(density(draws0), col="red")
curve(f, -3, 25, add=TRUE, col="black", n=1000)
(B0 <- coda::effectiveSize(draws0))   # Effective number of independent samples
B0 / time0['user.self']             # Effective number of independent samples per CPU second
counter0                           # Number of function evaluations

x_test = 0.7

integrate(f, -Inf, x_test)
mean(draws0 < x_test)
(err0 <- abs( integrate(f, -Inf, x_test)$value - mean(draws0 < x_test) ))

(acf0 <- acf(draws0)$acf[2,1,1])


# Transform
hist(draws1[-c(1:10000)], freq=FALSE, breaks=50)
lines(density(draws1), col="red")
curve(f, -3, 25, add=TRUE, col="black", n=1000)
curve(exp(pseudoLogPDF(x)), -3, 25, add=TRUE, col="blue", lty=2)
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


# Elliptical
hist(draws3[-c(1:10000)], freq=FALSE, breaks=50)
lines(density(draws3), col="red")
curve(f, -3, 25, add=TRUE, col="black", n=1000)
curve( dt((x-ellip_ctr)/ellip_sc, df=ellip_degf, log=TRUE) - log(ellip_sc), -3, 25, add=TRUE, col="blue", lty=2) # only makes sense if no transformation
(B3 <- coda::effectiveSize(draws3))  # Effective number of independent samples
B3 / time3['user.self']             # Effective number of independent samples per CPU second
counter3                            # Number of function evaluations

integrate(f, -Inf, x_test)
mean(draws3 < x_test)
(err3 <- abs( integrate(f, -Inf, x_test)$value - mean(draws3 < x_test) ))

(acf3 <- acf(draws3)$acf[2,1,1])





## compare
err0 # stepping out
err1 # transform
err2 # latent
err3 # t-elliptical

acf0
acf1
acf2
acf3

(B0 <- coda::effectiveSize(draws0))   # stepping out
(B1 <- coda::effectiveSize(draws1))   # transform
(B10 <- coda::effectiveSize(draws10))   # transform
(B2 <- coda::effectiveSize(draws2))   # latent
(B3 <- coda::effectiveSize(draws3))   # t-elliptical

B0 / time0['user.self'] # stepping out
B1 / time1['user.self'] # transform
B10 / time10['user.self'] # transform
B2 / time2['user.self'] # latent
B3 / time3['user.self'] # t-elliptical

counter0
counter1
counter2
counter3

