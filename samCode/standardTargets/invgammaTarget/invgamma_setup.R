## Set up for the Inverse gamma
## author: Sam Johnnson

filesToSource <- Sys.glob('../../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)

library(magrittr)
library(tidyverse)
library(LaplacesDemon, include.only = 'Thin') # thin function, KLD, JSD
library(cucumber)
library(utils)
library(rootSolve)


reps <- 100

fexp <- function(f, x) exp(f(x))
samples <- rep(5e4,reps)
saveInd <- TRUE


# 1
################### TEST ########################
## dinvgamma(x, shape = 2, scale = 1) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  dinvgamma(x, shape = 2, scale = 1, log = TRUE)
}


# making a grid to calculate KL divergence
# grid <- seq(from = 0,
#             to = qinvgamma((1-.99999)/2, shape = 2, scale = 1, lower.tail = FALSE),
#             length.out = 5000)
# 
# py <- exp(lf(grid))

truth = list(d = function(x) {dinvgamma(x, shape = 2, scale = 1)},
             ld = function(x) {dinvgamma(x, shape = 2, scale = 1,  log=TRUE)},
             # dld = function(x) {1.5/x - 1.0},
             q = function(u) {qinvgamma(u, shape = 2, scale = 1)},
             t = "invgamma(2,1)",
             lb = 0,
             ub = Inf)
truth$dld = function(x) numDeriv::grad(truth$ld, x=x)

