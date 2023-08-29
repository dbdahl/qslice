attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)

fixedBeta <- c(-0.02742211, 0.25728077, -0.23017904, 0.06894651, -0.56840819,
               0.23139456, 0.02402770,  0.19669150, 0.07542407, -0.05148203)
fixedPsi <- 5.264618

# beta ~ mvn(beta_0, g / psi * inv(XtX)), using the covariance parametrization
yty = t(y) %*% y
Xt <- t(X)
XtX <- Xt %*% X
inv_XtX <- solve(XtX)
beta_mle <- solve(Xt %*% X, Xt %*% y)
beta_0 <- rep(0, ncol(X))
beta <- fixedBeta

# psi ~ gamma(a, b)
a_0 <- 5
b_0 <- 1
psi <- a_0 / b_0

# g ~ uniform(0, g_max)
g_max <- 2 * ncol(X)^2
g <- g_max / 2

f <- function(g) {g^(-n/2) * exp(-0.5 * (psi/g) * t(beta-beta_0)%*%XtX%*%(beta-beta_0))}
c <- 0.5 * psi * t(beta-beta_0)%*%XtX%*%(beta-beta_0)
mode <- (2*c)/(n)
sigma <- function(g) {
  x1 <- (2*c)/(g^3)
  x2 <- (n)/(2*g^2)
  x1x2 <- x1 - x2
  sqrt(1/x1x2)
}

curve(dnorm(x,mean = mode, sd =  sigma(mode))/0.45, 0,20, col = 'blue')
curve(f(x)/5e-17,0,20, add = TRUE)
abline(v = mode, col = 'red')
