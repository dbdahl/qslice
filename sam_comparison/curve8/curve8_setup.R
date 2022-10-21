source('../setup.R')
curve_num <- 8

# 8
################### Curve 8 #############################
## (x^3 + x^2 + 100*sin(x))(-100*cos(5) + (3575)/(12)) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  if_else(x < 0, log(0),
          if_else(x <= 5, log((x^3 + x^2 + 100*sin(x))/(-100*cos(5) + (3575)/(12))),
                  log(0)))
}

cdf <- function(x) {
  probs <- vector(length = length(x))
  for (i in 1:length(x)) {
    probs[i] <- if_else(x <0, 0,
                        if_else(x <= 5, ((12)/(-1200*cos(5) + 3575))*((x^4)/(4) + (x^3)/(3) + 100*(-cos(x) + 1)), 1))
  }
  return(probs)
}

# plot of the log of the density function below
pdf(file = "../../images_slice_sampler_comp/curve8.pdf")
curve(fexp(f = lf, x = x), xlim = c(0,5), ylim = c(0, 0.4))
dev.off()

grid <- seq(from = 0,
            to = 5,
            length.out = 1000)

py <- exp(lf(grid))

xlim_range <- c(0,3)
ylim_range <- c(0, 0.4)

#### Tuning Parameters ####
## starting point ##
x <- c(0.5, 1, 4.5)
## stepping out metrics to input ##
w <- c(2, 5, 10)

## latent slice sampling metric to input ##
s <- c(3, 5, 10)
rate <- c(2)

## gess slice sampling metrics to input ##
mu <- c(3, 5, 10)
sigma <- c(3)
df <- c(3)

## transform tuning parameters ##
log_pdf <- c(function(x) dnorm(x, mean = 5, sd = 5, log = TRUE),
             function(x) dnorm(x, mean = 4, sd = 10, log = TRUE))

inv_cdf <- c(function(u) qnorm(u, mean = 5, sd = 5),
             function(u) qnorm(u, mean = 4, sd = 10))

find_grid <- list(seq(from = qnorm((1-.99999)/2, mean = 5, sd = 5),
                      to = qnorm((1-.99999)/2, mean = 5, sd = 5, lower.tail = FALSE), length.out = 1000),
                  seq(from = qnorm((1-.99999)/2, mean = 4, sd = 10),
                      to = qnorm((1-.99999)/2, mean = 4, sd = 10, lower.tail = FALSE), length.out = 1000))


## random walk tuning parameters ##
c <- c(1)
