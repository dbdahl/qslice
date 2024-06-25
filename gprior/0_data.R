data(mtcars)

dat <- list()

attach(mtcars)
dat$y <- scale(mpg)
dat$X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
detach(mtcars)

(dat$n <- length(dat$y))
(dat$p <- ncol(dat$X))

dat$XtX <- crossprod(dat$X)
(dat$logdet_XtX <- determinant(dat$XtX, log = TRUE)$modulus)
dat$chol_XtX <- chol(dat$XtX)
dat$inv_chol_XtX <- backsolve(dat$chol_XtX, diag(dat$p))
(dat$beta_mle <- solve(dat$XtX, crossprod(dat$X, dat$y)) |> drop())

dat
rm(mtcars)
