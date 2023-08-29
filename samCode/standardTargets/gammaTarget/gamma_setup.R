## Set up for the cauch transform 
## author: Sam Johnnson

# setwd("~/cucumber/sam_comparison")
filesToSource <- Sys.glob('../../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)

library(magrittr)
library(tidyverse)
library(LaplacesDemon) # thin function, KLD, JSD
library(cucumber)
library(utils)
library(rootSolve)

reps <- 100

fexp <- function(f, x) exp(f(x))
samples <- rep(5e4,reps)
saveInd <- TRUE


# 1
################### TEST ########################
## dgamma(x,2.5,1) ##

# curve 1 is bimodal with one mode being much larger
lf <- function(x) {
  dgamma(x,2.5,1,log = TRUE)
}

truth = list(d = function(x) {dgamma(x, 2.5, 1)},
             ld = function(x) {dgamma(x, 2.5, 1,  log=TRUE)},
             # dld = function(x) {1.5/x - 1.0},
             q = function(u) {qgamma(u, 2.5, 1)},
             t = "gamma(2.5,1)",
             lb = 0,
             ub = Inf)
truth$dld = function(x) numDeriv::grad(truth$ld, x=x)

xlim_range <- c(-4, 15)
ylim_range <- c(0, 0.42)


