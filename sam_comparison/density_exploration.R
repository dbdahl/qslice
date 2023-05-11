## This script explores the density function ##
## author: Sam Johnson ##

# creating functions to use
# plots log scale on the normal
fexp <- function(f, x) exp(f(x))
# draws from curve 1
curve1Sampler <- function() {
  unif <- runif(1,0,1)
  ifelse(unif < 1/5, rnorm(1, mean = 0, sd = 0.5), rnorm(1, mean = 6, sd = 2))
}

cdf <- function(x){
  0.2*pnorm(x, sd = 0.5) + 0.8 * pnorm(x, mean = 6, sd =2)
}
# curve 1
lf <- function(x) {
  log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 6, sd = 2))
}
# curve 1 cdf
# cdf <- function(x, mean2 = 6, sd2 = 2) {
#   probs <- vector(length = length(x))
#   for (i in 1:length(x)) {
#     probs[i] <-
#       integrate(function(x) {
#         exp(log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = mean2, sd = sd2)))
#       },
#       lower = -Inf, upper = x[i])$value
#   }
#   return(probs)
# }

# taking 1e6 draws from curve1
curve1Draws <- replicate(1e6, curve1Sampler())
# plotting density curve
# bw options: nrd0, nrd, ucv, bcv, SJ
# bwOptions <- c('nrd0','nrd','ucv','bcv','SJ')
# par(mfrow = c(2,3))
# for(i in 1:length(bwOptions)) {
#   plot(density(curve1Draws, bw = bwOptions[i]), main = bwOptions[i]) 
#   curve(fexp(lf,x), col = 'red', add = TRUE, type = 'l')
# }

ks.test(curve1Draws, cdf)

curve(fexp(lf,x), col = 'red', type = 'l', from =-5, to = 15)
lines(density(curve1Draws), col = 'dodgerblue')

curve4 <- function(x) {
  ifelse(x < 0, -Inf, dt(x, df = 5.0, log = TRUE) + log(2.0))
}

curve(fexp(curve4, x), from = -1, to = 5)
