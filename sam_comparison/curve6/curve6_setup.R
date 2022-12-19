source('../setup.R')
curve_num <- 6

# 6
################### Curve 6 ####################
## dbeta(x, shape1=0.2, shape2=0.8, log=TRUE) ##

# u shaped distribution
lf <- function(x) {
  dbeta(x,
        shape1 = 0.2,
        shape2 = 0.8,
        log = TRUE)
}
# plot of the log of the density function below
pdf(file = "../../images_slice_sampler_comp/curve6.pdf")
curve(fexp(f = lf, x = x), xlim = c(0, 1), ylim = c(0, 10))
dev.off()

grid <- seq(from = 0.00001,
            to = 0.99999,
            length.out = 1000)

# value for Kullback-Leibler Divergence
py <- exp(lf(grid))

xlim_range <- c(0, 1)
ylim_range <- c(0, 5 + 0.10)

#### Tuning Parameters ####
## starting point ##
x <- c(0.2)#c(0.2, 0.5, 0.8)

## stepping out metrics to input ##
w <- c(0.4, 1, 2, 5, 10)#c(0.01, 1, 2, 4, 10)

## latent slice sampling metric to input ##
s <- c(3)#c(3, 5, 10)#c(0.01, 1, 2, 10)
rate <- c(0.1, 0.5, 2, 5, 10)#c(0.5, 1, 1.5, 2, 2.5, 3)

## gess slice sampling metrics to input ##
mu <- c(0.5,0.8)#c(1,2,3,4.5,6,7)
sigma <- c(0.25, 0.5, 0.75)#c(2,3,4,5,6,8)
df <- c(3, 5, 10, 15)#c(1,4,16,16^2,16^4)

## transform tuning parameter ##
log_pdf <- c(function(x) dbeta(x, shape1 = 0.3, shape2 = 0.7, log = TRUE),
             function(x) dunif(x, min = 0, max = 1, log = TRUE),
             function(x) dbeta(x, shape1 = 0.5, shape2 = 0.5, log = TRUE),
             function(x) dbeta(x, shape1 = 5, shape2 = 5, log = TRUE),
             function(x) dbeta(x, shape1 = 0.1, shape2 = 0.9, log = TRUE),
             function(x) dbeta(x, shape1 = 4, shape2 = 5, log = TRUE),
             function(x) dbeta(x, shape1 = 2, shape2 = 1, log = TRUE),
             function(x) dbeta(x, shape1 = 1, shape2 = 1, log = TRUE),
             function(x) dbeta(x, shape1 = 9, shape2 = 1, log = TRUE),
             function(x) dbeta(x, shape1 = 2, shape2 = 2, log = TRUE)
             )

inv_cdf <- c(function(u) qbeta(u, shape1 = 0.3, shape2 = 0.7),
             function(u) qunif(u, min = 0, max = 1),
             function(u) qbeta(u, shape1 = 0.5, shape2 = 0.5),
             function(u) qbeta(u, shape1 = 5, shape2 = 5),
             function(u) qbeta(u, shape1 = 0.1, shape2 = 0.9),
             function(u) qbeta(u, shape1 = 4, shape2 = 5),
             function(u) qbeta(u, shape1 = 2, shape2 = 1),
             function(u) qbeta(u, shape1 = 1, shape2 = 1),
             function(u) qbeta(u, shape1 = 9, shape2 = 1),
             function(u) qbeta(u, shape1 = 2, shape2 = 2)
             )

find_grid <- list(seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000),
                  seq(from = 0.00001,
                      to = 0.99999, length.out = 1000)
                  )

## rand walk tuning parameter ##
c <- c(0.1, 0.5, 1)
