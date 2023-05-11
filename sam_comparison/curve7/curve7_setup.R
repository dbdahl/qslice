source('../setup.R')
curve_num <- 7

# 7
############## Curve 7 ######################
## dnorm(x, mean = 20, sd = 5, log = TRUE) ##

# normal prior normal likelihood function
lf <- function(x) {
  dnorm(x, mean = 20, sd = 5, log = TRUE)
}
# plot of the density !! can use for elliptical slice sampler !!
pdf(file = "../../images_slice_sampler_comp/curve7.pdf")
curve(fexp(f = lf, x = x), xlim = c(0, 40))
dev.off()

grid <- seq(from = qnorm((1-.9999)/2, mean = 20, sd = 5),
            to = qnorm((1-.9999)/2, mean = 20, sd = 5, lower.tail = FALSE),
            length.out = 5000)

# value for Kullback-Leibler Divergence
py <- exp(lf(grid))

xlim_range <- c(0, 40)
ylim_range <- c(0, 0.08 + 0.010)

#### Tuning Parameters ####
## starting point ##
x <- c(1)#c(0, 1, 5)

## stepping out metrics to input ##
w <- c(0.5, 1, 2, 5, 10, 20, 30)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.1, 0.5, 2, 5, 10)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(20, 15, 25)#c(1,2,3,4.5,6,7)
sigma <- c(5, 3, 8)#c(2,3,4,5,6,8)
df <- c(3)#c(1,4,16,16^2,16^4)

## transform tuning parameters ##
laplace_approximation <- laplace_approx(lf, init = 1)

log_pdf <- c(function(x) dnorm(x, mean = 20, sd = 6, log = TRUE),
             function(x) dnorm(x, mean = 15, sd = 8, log = TRUE),
             laplace_approximation$log_pdf,
             function(x) dnorm(x, mean = 20, sd = 5, log = TRUE),
             function(x) dnorm(x, mean = 20, sd = 4, log = TRUE),
             function(x) dnorm(x, mean = 15, sd = 10, log = TRUE),
             function(x) dnorm(x, mean = 20, sd = 10, log = TRUE),
             function(x) dnorm(x, mean = 15, sd = 15, log = TRUE),
             function(x) dnorm(x, mean = 10, sd = 20, log = TRUE),
             function(x) dnorm(x, mean = 0, sd = 30, log = TRUE)
             )

inv_cdf <- c(function(u) qnorm(u, mean = 20, sd = 6),
             function(u) qnorm(u, mean = 15, sd = 8),
             laplace_approximation$inv_cdf,
             function(u) qnorm(u, mean = 20, sd = 5),
             function(u) qnorm(u, mean = 20, sd = 4),
             function(u) qnorm(u, mean = 15, sd = 10),
             function(u) qnorm(u, mean = 20, sd = 10),
             function(u) qnorm(u, mean = 15, sd = 15),
             function(u) qnorm(u, mean = 10, sd = 20),
             function(u) qnorm(u, mean = 0, sd = 30)
             )

## rand walk tuning parameters ##
c <- c(1, 2, 5, 10)
