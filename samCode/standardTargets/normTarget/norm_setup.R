## Set up for the cauch transform 
## author: Sam Johnnson

filesToSource <- Sys.glob('../../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)

library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)
library(rootSolve)
library(parallel)


reps <- 100

fexp <- function(f, x) exp(f(x))
samples <- rep(5e4,reps)
saveInd <- TRUE


# 1
################### TEST ########################
## dnorm(x,0,1) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  dnorm(x,0,1,log = TRUE)
}

truth = list(d = function(x) {dnorm(x)},
             ld = function(x) {dnorm(x, log=TRUE)},
             # dld = function(x) {-x},
             q = function(u) {qnorm(u)},
             lb = -Inf, ub = Inf,
             t = "normal(0,1)")
truth$dld = function(x) numDeriv::grad(truth$ld, x=x)

xlim_range <- c(-3, 3)
ylim_range <- c(0, 0.42)
