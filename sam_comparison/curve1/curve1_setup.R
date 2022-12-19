# setwd("~/cucumber/sam_comparison")
source('../setup.R')
curve_num <- 1

# 1
################### Curve 1 ########################
## 0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=6,sd=2) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 6, sd = 2))
}

cdf <- function(x, mean2 = 6, sd2 = 2) {
  probs <- vector(length = length(x))
  for (i in 1:length(x)) {
    probs[i] <-
      integrate(function(x) {
        exp(log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = mean2, sd = sd2)))
      },
      lower = -Inf, upper = x[i])$value
  }
  return(probs)
}

# plot of the log of the density function below
pdf(file = "../../images_slice_sampler_comp/curve1.pdf")
curve(fexp(f = lf, x = x), xlim = c(-4, 15), ylim = c(0, .17))
dev.off()

grid <- seq(from = qnorm((1-.99999)/2, mean = 0, sd = 0.5),
            to = qnorm((1-.99999)/2, mean = 6, sd = 2, lower.tail = FALSE),
            length.out = 1000)

py <- exp(lf(grid))

xlim_range <- c(-4, 15)
ylim_range <- c(0, 0.17 + 0.10)

#### Tuning Parameters ####
## starting point ##
x <- c(1)#c(0, 1, 5)

## stepping out metrics to input ##
w <- c(0.5, 1, 2, 5, 10, 15, 20)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.1, 0.5, 1, 2, 5)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(3, 5, 10)#c(1,2,3,4.5,6,7)
sigma <- c(2, 3, 6)#c(2,3,4,5,6,8)
df <- c(3, 5, 30)#c(1,4,16,16^2,16^4)

## transform tuning parameters ##
log_pdf <- c(function(x) dnorm(x, mean = 5, sd = 5, log = TRUE),
             function(x) dnorm(x, mean = 4, sd = 10, log = TRUE),
             function(x) dnorm(x, mean = 0, sd = 6, log = TRUE),
             function(x) dnorm(x, mean = 6, sd = 6, log = TRUE),
             function(x) dnorm(x, mean = 3, sd = 6, log = TRUE),
             function(x) dnorm(x, mean = 1, sd = 8, log = TRUE),
             function(x) dnorm(x, mean = 6, sd = 11, log = TRUE),
             function(x) dnorm(x, mean = 6, sd = 3, log = TRUE),
             function(x) dnorm(x, mean = 2, sd = 11, log = TRUE),
             function(x) dnorm(x, mean = 7, sd = 8, log = TRUE)
             )

inv_cdf <- c(function(u) qnorm(u, mean = 5, sd = 5),
             function(u) qnorm(u, mean = 4, sd = 10),
             function(u) qnorm(u, mean = 0, sd = 6),
             function(u) qnorm(u, mean = 6, sd = 6),
             function(u) qnorm(u, mean = 3, sd = 6),
             function(u) qnorm(u, mean = 1, sd = 9),
             function(u) qnorm(u, mean = 6, sd = 11),
             function(u) qnorm(u, mean = 6, sd = 3),
             function(u) qnorm(u, mean = 2, sd = 11),
             function(u) qnorm(u, mean = 7, sd = 8)
             )

find_grid <- list(seq(from = qnorm((1-.99999)/2, mean = 5, sd = 5),
                      to = qnorm((1-.99999)/2, mean = 5, sd = 5, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 4, sd = 10),
                      to = qnorm((1-.99999)/2, mean = 4, sd = 10, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 6),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 6, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 6, sd = 6),
                      to = qnorm((1-.99999)/2, mean = 6, sd = 6, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 3, sd = 6),
                      to = qnorm((1-.99999)/2, mean = 3, sd = 6, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 1, sd = 9),
                      to = qnorm((1-.99999)/2, mean = 1, sd = 9, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 6, sd = 11),
                      to = qnorm((1-.99999)/2, mean = 6, sd = 11, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 6, sd = 3),
                      to = qnorm((1-.99999)/2, mean = 6, sd = 3, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 2, sd = 11),
                      to = qnorm((1-.99999)/2, mean = 2, sd = 11, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 7, sd = 8),
                      to = qnorm((1-.99999)/2, mean = 7, sd = 8, lower.tail = FALSE), length.out = 1000)
                  )

## random walk tuning parameters ##
c <- c(1,2,5,10)
