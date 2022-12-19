# setwd("~/cucumber/sam_comparison")
source('../setup.R')
curve_num <- 2

# 2
################## Curve2 ###########################
## 0.2*dnorm(x,sd=0.5) + 0.8*dnorm(x,mean=20,sd=1) ##

# curve 2 is bimodal
lf <- function(x) {
  log(0.2 * dnorm(x, sd = 0.5) + 0.8 * dnorm(x, mean = 20, sd = 1))
}

cdf <- function(x, mean2 = 20, sd2 = 1) {
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
pdf(file = "../../images_slice_sampler_comp/curve2.pdf")
curve(fexp(f = lf, x = x), xlim = c(-3, 25), ylim = c(0, .32))
dev.off()

grid <- seq(from = qnorm((1-.99999)/2, mean = 0, sd = 0.5),
            to = qnorm((1-.99999)/2, mean = 20, sd = 1, lower.tail = FALSE),
            length.out = 1000)

py <- exp(lf(grid))

xlim_range <- c(-3, 25)
ylim_range <- c(0, 0.32 + 0.10)

#### Tuning Parameters ####
## starting point ##
x <- c(0)#c(0, 5, 10)

## stepping out metrics to input ##
w <- c(1, 2, 5, 10, 20)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.1, 0.5, 2, 5)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(3, 5, 10)#c(1,2,3,4.5,6,7)
sigma <- c(3, 5, 10)#c(2,3,4,5,6,8)
df <- c(3, 5, 10)#c(1,4,16,16^2,16^4)

## transform slice sampling tuning parameters ##

log_pdf <- c(
  function(x) dnorm(x, mean = 10, sd = 10, log = TRUE),
  function(x) dnorm(x, mean = 0, sd = 20, log = TRUE),
  function(x) dnorm(x, mean = 10, sd = 5, log = TRUE),
  function(x) dnorm(x, mean = 10, sd = 20, log = TRUE),
  function(x) dnorm(x, mean = 10, sd = 30, log = TRUE),
  function(x) dnorm(x, mean = 20, sd = 20, log = TRUE),
  function(x) dnorm(x, mean = 15, sd = 20, log = TRUE),
  function(x) dnorm(x, mean = 20, sd = 30, log = TRUE),
  function(x) dnorm(x, mean = 0, sd = 30, log = TRUE),
  function(x) dnorm(x, mean = 7, sd = 30, log = TRUE)
)

inv_cdf <- c(function(u) qnorm(u, mean = 10, sd = 10),
             function(u) qnorm(u, mean = 0, sd = 20),
             function(u) qnorm(u, mean = 10, sd = 5),
             function(u) qnorm(u, mean = 10, sd = 20),
             function(u) qnorm(u, mean = 10, sd = 30),
             function(u) qnorm(u, mean = 20, sd = 20),
             function(u) qnorm(u, mean = 15, sd = 20),
             function(u) qnorm(u, mean = 20, sd = 30),
             function(u) qnorm(u, mean = 0, sd = 30),
             function(u) qnorm(u, mean = 7, sd = 30)
)

find_grid <- list(seq(from = qnorm((1-.99999)/2, mean = 10, sd = 10),
                      to = qnorm((1-.99999)/2, mean = 10, sd = 10, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 20),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 20, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 10, sd = 5),
                      to = qnorm((1-.99999)/2, mean = 10, sd = 5, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 10, sd = 20),
                      to = qnorm((1-.99999)/2, mean = 10, sd = 20, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 10, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 10, sd = 30, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 20, sd = 20),
                      to = qnorm((1-.99999)/2, mean = 20, sd = 20, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 15, sd = 20),
                      to = qnorm((1-.99999)/2, mean = 15, sd = 20, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 20, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 20, sd = 30, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 0, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 0, sd = 30, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 7, sd = 30),
                      to = qnorm((1-.99999)/2, mean = 7, sd = 30, lower.tail = FALSE), length.out = 1000)
                  )

## random walk tuning parameters ##
c <- c(1, 2, 5, 10)