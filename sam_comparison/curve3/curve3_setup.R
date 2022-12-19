source('../setup.R')
curve_num <- 3

# 3
############ Curve 3 ###############
## dgamma(x, shape=2.5, log=TRUE) ##

# curve 3 in unimodal skewed right
lf <- function(x) {
  dgamma(x,
         shape = 2.5,
         rate = 1,
         log = TRUE)
}
# plot of the log of the density function below
pdf(file = "../../images_slice_sampler_comp/curve3.pdf")
curve(fexp(f = lf, x = x), xlim = c(0, 10), ylim = c(0, .31))
dev.off()

grid <- seq(from = 0,
            to = qgamma(.99999, shape = 2.5, rate = 1),
            length.out = 1000)

# value for Kullback-Leibler Divergence
py <- exp(lf(grid))

xlim_range <- c(0, 10)
ylim_range <- c(0, 0.31 + 0.10)

#### Tuning Parameters ####
## starting point ##
x <- c(1)#c(0.5, 1, 5)

## stepping out metrics to input ##
w <- c(1, 2, 3, 5, 10, 20)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.1, 0.5, 2, 5)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(3, 5, 10)#c(1,2,3,4.5,6,7)
sigma <- c(1, 3, 5)#c(2,3,4,5,6,8)
df <- c(3, 5, 10)#c(1,4,16,16^2,16^4)

## transformation tuning parameters ##
laplace_approximation <- laplace_approx(lf, init = 1)

log_pdf <- c(function(x) dgamma(x, shape = 2.5, rate = 1, log = TRUE),
             function(x) dnorm(x, mean = 2, sd = 4, log = TRUE),
             laplace_approximation$log_pdf,
             function(x) dnorm(x, mean = 4, sd = 4, log = TRUE),
             function(x) dgamma(x, shape = 2.5, rate = 2, log = TRUE),
             function(x) dgamma(x, shape = 2.5, rate = 3, log = TRUE),
             function(x) dgamma(x, shape = 10, rate = 4, log = TRUE),
             function(x) dgamma(x, shape = 3, rate = 1, log = TRUE),
             function(x) dgamma(x, shape = 3.5, rate = 1, log = TRUE),
             function(x) dgamma(x, shape = 4, rate = 1, log = TRUE)
)

inv_cdf <- c(function(u) qgamma(u, shape = 2.5, rate = 1),
             function(u) qnorm(u, mean = 2, sd = 4),
             laplace_approximation$inv_cdf,
             function(u) qnorm(u, mean = 4, sd = 4),
             function(u) qgamma(u, shape = 2.5, rate = 2),
             function(u) qgamma(u, shape = 2.5, rate = 3),
             function(u) qgamma(u, shape = 10, rate = 4),
             function(u) qgamma(u, shape = 3, rate = 1),
             function(u) qgamma(u, shape = 3.5, rate = 1),
             function(u) qgamma(u, shape = 4, rate = 1)
)

find_grid <- list(seq(from = 0, to = qgamma(0.99999, shape = 2.5, rate = 1), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 2, sd = 4),
                      to = qnorm((1-.99999)/2, mean = 2, sd = 4, lower.tail = FALSE), length.out = 1000),
                  seq(from = laplace_approximation$inv_cdf((1-0.99999)/2),
                      to = laplace_approximation$inv_cdf((1-0.99999)/2, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 4, sd = 4),
                      to = qnorm((1-.99999)/2, mean = 4, sd = 4, lower.tail = FALSE), length.out = 1000),
                  seq(from = 0, to = qgamma(0.99999, shape = 2.5, rate = 2), length.out = 1000),
                  seq(from = 0, to = qgamma(0.99999, shape = 2.5, rate = 3), length.out = 1000),
                  seq(from = 0, to = qgamma(0.99999, shape = 10, rate = 4), length.out = 1000),
                  seq(from = 0, to = qgamma(0.99999, shape = 3, rate = 1), length.out = 1000),
                  seq(from = 0, to = qgamma(0.99999, shape = 3.5, rate = 1), length.out = 1000),
                  seq(from = 0, to = qgamma(0.99999, shape = 4, rate = 1), length.out = 1000)
                  )

## random walk tuning parameters ##
c <- c(1, 2, 5, 10)
