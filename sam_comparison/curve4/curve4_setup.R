source('../setup.R')
curve_num <- 4

# 4
######################### Curve 4 ############################
## ifelse( x < 0, -Inf, dt(x, df=3.0, log=TRUE) + log(2.0)) ##

# curve 4 is right skewed with support for x strictly positive
lf <- function(x) {
  ifelse(x < 0, -Inf, dt(x, df = 3.0, log = TRUE) + log(2.0))
}

cdf <- function(x, df1 = 3.0) {
  probs <- vector(length = length(x))
  for (i in 1:length(x)) {
    probs[i] <-
      integrate(function(x) {
        exp(ifelse(x < 0, -Inf, dt(x, df = df1, log = TRUE) + log(2.0)))
      },
      lower = 0, upper = x[i])$value
  }
  return(probs)
}

# plot of the log of the density function below
pdf(file = "../../images_slice_sampler_comp/curve4.pdf")
curve(fexp(f = lf, x = x), xlim = c(0, 10), ylim = c(0, .8))
dev.off()

grid <- seq(from = 0,
            to = qt(.99999, df = 3),
            length.out = 1000)

# value for Kullback-Leibler Divergence
py <- exp(lf(grid))

xlim_range <- c(0, 10)
ylim_range <- c(0, 0.8 + 0.10)

#### Tuning Parameters ####
## starting point ##
x <- c(1)#c(0, 1, 5)

## stepping out metrics to input ##
w <- c(0.5, 1, 2, 5, 10, 20)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.1, 0.5, 1, 2, 5)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(3, 5, 10)#c(1,2,3,4.5,6,7)
sigma <- c(1, 3, 5)#c(2,3,4,5,6,8)
df <- c(3, 5, 10)#c(1,4,16,16^2,16^4)

## transform tuning parameters ##
log_pdf <- c(function(x) dt(x, df = 3, log = TRUE),
             function(x) dnorm(x, mean = 0, sd = 15, log = TRUE),
             function(x) dt(x, df = 2, log = TRUE),
             function(x) dt(x, df = 1, log = TRUE),
             function(x) dnorm(x, mean = 2, sd = 15, log = TRUE),
             function(x) dnorm(x, mean = 0, sd = 25, log = TRUE),
             function(x) dnorm(x, mean = 2, sd = 25, log = TRUE),
             function(x) dnorm(x, mean = 0, sd = 30, log = TRUE),
             function(x) dnorm(x, mean = 2, sd = 30, log = TRUE),
             function(x) dnorm(x, mean = 4, sd = 30, log = TRUE)
             )

inv_cdf <- c(function(u) qt(u, df = 3),
             function(u) qnorm(u, mean = 0, sd = 15),
             function(u) qt(u, df = 2),
             function(u) qt(u, df = 1),
             function(u) qnorm(u, mean = 2, sd = 15),
             function(u) qnorm(u, mean = 0, sd = 25),
             function(u) qnorm(u, mean = 2, sd = 25),
             function(u) qnorm(u, mean = 0, sd = 30),
             function(u) qnorm(u, mean = 2, sd = 30),
             function(u) qnorm(u, mean = 4, sd = 30)
             )

find_grid <- list(seq(from = qt((1-.99999)/2, df = 3),
                      to = qt((1-.99999)/2, df = 3, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 3),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 3, lower.tail = FALSE), length.out = 1000),
                  seq(from = qt((1-.99999)/2, df = 2),
                      to = qt((1-.99999)/2, df = 2, lower.tail = FALSE), length.out = 1000),
                  seq(from = qt((1-.99999)/2, df = 1),
                      to = qt((1-.99999)/2, df = 1, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 2, sd = 15),
                      to = qnorm((1-.99999)/2, mean = 2, sd = 15, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 25),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 25, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 2, sd = 25),
                      to = qnorm((1-.99999)/2, mean = 2, sd = 25, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 30, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 2, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 2, sd = 30, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 4, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 4, sd = 30, lower.tail = FALSE), length.out = 1000)
                  )

## random walk tuning parameters ##
c <- c(1, 2, 3, 5, 10)
