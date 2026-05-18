
dat <- list()

if (data_use == "mtcars") {

  data(mtcars)

  attach(mtcars)
  dat$y <- scale(mpg)
  dat$X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
  detach(mtcars)

  rm(mtcars)

  dat$XtX <- crossprod(dat$X)
  (dat$beta_hat <- solve(dat$XtX, crossprod(dat$X, dat$y)) |> drop())

} else if (grep("db", data_use)) {

  library("lars")
  data(diabetes)
  str(diabetes)
  # ?diabetes

  if (data_use == "db40") {
    set.seed(1)
    indx_obs <- sort(sample.int(length(diabetes$y), replace = FALSE, size = 40))
  } else {
    indx_obs <- 1:length(diabetes$y)
  }

  dat$y <- scale(diabetes$y[indx_obs])
  dat$X <- scale(diabetes$x2[indx_obs,])

  rm(diabetes)

  dat$XtX <- crossprod(dat$X)
  library("glmnet")
  lso <- cv.glmnet(dat$X, dat$y, alpha = 1, family = "gaussian")
  dat$beta_hat <- as.numeric(coef(lso)[-1])
  str(dat$beta_hat)

} else if (data_use == "riboflavin") {

  library("hdi")
  data(riboflavin)
  str(riboflavin)
  # hist(riboflavin$y)
  # hist(riboflavin$x[,1])
  # ?riboflavin

  dat$y <- scale(riboflavin$y)
  dat$X <- scale(riboflavin$x)

  rm(riboflavin)

  dat$XtX <- crossprod(dat$X)

  library("glmnet")
  lso <- cv.glmnet(dat$X, dat$y, alpha = 1, family = "gaussian")
  dat$beta_hat <- as.numeric(coef(lso)[-1])
  str(dat$beta_hat)

}


## Common to all data sets

(dat$n <- length(dat$y))
(dat$p <- ncol(dat$X))

dat
