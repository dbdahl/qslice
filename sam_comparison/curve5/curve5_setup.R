source('../setup.R')

curve_num <- 5

# 5
########## Curve 5 ##########
## dt(x, df=5.0, log=TRUE) ##

# t distribution with a steep peak
lf <- function(x) {
  dt(x, df = 5.0, log = TRUE)
}
# plot of the log of the density function below
pdf(file = "../../images_slice_sampler_comp/curve5.pdf")
curve(fexp(f = lf, x = x), xlim = c(-10, 10), ylim = c(0, .42))
dev.off()

grid <- seq(from = qt((1-.99999)/2, df = 5),
            to = qt((1-.99999)/2, df = 5, lower.tail = FALSE),
            length.out = 1000)

# value for Kullback-Leibler Divergence
py <- exp(lf(grid))


xlim_range <- c(-10, 10)
ylim_range <- c(0, 0.42 + 0.10)

#### Tuning Parameters ####
## starting point ##
x <- c(0, 1, 5)

## stepping out metrics to input ##
w <- c(2, 5, 10)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(2)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(0, 5, 10)#c(1,2,3,4.5,6,7)
sigma <- c(3)#c(2,3,4,5,6,8)
df <- c(1,3)#c(1,4,16,16^2,16^4)

## tranform tuning parameter ##
laplace_approximation <- laplace_approx(lf, init = 1)

log_pdf <- c(function(x) dt(x, df = 2, log = TRUE),
             function(x) dnorm(x, mean = 0, sd = 8, log = TRUE),
             laplace_approximation$log_pdf)

inv_cdf <- c(function(u) qt(u, df = 2),
             function(u) qnorm(u, mean = 0, sd = 8),
             laplace_approximation$inv_cdf)

find_grid <- list(seq(from = qt((1-.99999)/2, df = 2),
                      to = qt((1-.99999)/2, df = 2, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 8),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 8, lower.tail = FALSE), length.out = 1000),
                  seq(from = laplace_approximation$inv_cdf((1-0.99999)/2),
                      to = laplace_approximation$inv_cdf((1-0.99999)/2, lower.tail = FALSE), length.out = 1000))

## random walk tuning parameter ##
c <- c(5,8)
