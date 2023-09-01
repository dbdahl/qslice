## gprior setup function
## author: Sam Johnson

library(mvtnorm)
library(cucumber)
library(parallel)
library(mcmc)
library(coda)

filesToSource <- Sys.glob('../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)

nChains <- 5
nSamples <- 55000
nThin <- 2
nBurnin <- 5000

data(mtcars)

attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)

# beta ~ mvn(beta_0, g / psi * inv(XtX)), using the covariance parametrization
yty = t(y) %*% y
Xt <- t(X)
XtX <- Xt %*% X
inv_XtX <- solve(XtX)
beta_mle <- solve(Xt %*% X, Xt %*% y)
beta_0 <- rep(0, ncol(X))
beta <- beta_0

# psi ~ gamma(a, b)
a_0 <- 5
b_0 <- 1
psi <- a_0 / b_0

# g ~ uniform(0, g_max)
g_max <- 2 * ncol(X)^2
g <- g_max / 2


# the samplers for each method
# stepping out
steppingOutList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                        method = 'SteppingOut', w = 30)
# latent
latentList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                   method = 'Latent', s = 5, rate = 0.001)
# gess
gessList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                 method = 'GESS', mu = 40, sigma = 30, df = 3)
# randwalk
randWalkList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                     method = 'RandWalk', c = 30)
# independence
independenceList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                         method = 'Independence', pseudoType = 'OptimSamplesAUC',
                         optimType = 'samples', approxSamples = NULL, everyIter = FALSE)
# transform
# auto
autoList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                 method = 'Transform', pseudoType = 'Auto',
                 approxSamples = NULL, everyIter = FALSE)

# laplace
laplaceList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                    method = 'Transform', pseudoType = 'Laplace', lapalceFunc = NULL,
                    everyIter = 1)
# optim samples
optimSamplesList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                         method = 'Transform', pseudoType = 'OptimSamples',
                         optimType = 'samples', approxSamples = NULL, everyIter = FALSE)

# optim samples
optimSamplesAUCList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                         method = 'Transform', pseudoType = 'OptimSamplesAUC',
                         optimType = 'samples', approxSamples = NULL, everyIter = FALSE)

# optim
optimList <- list(Nsamples = nSamples, Nburnin = nBurnin, Nthin = nThin,
                  method = 'Transform', pseudoType = 'Optim', optimType = 'function',
                  everyIter = FALSE)

