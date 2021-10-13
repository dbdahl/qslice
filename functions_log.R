stunif <- function() runif(1, 0, 1)

# Figure 3 of Neal (2003)
slice_sampler_using_stepping_out <- function(x, lf, w=1, max_number_of_steps=Inf) {
  # Step 1
  if ( is.list(x) ) {
    lfx <- x$lfx
    x <- x$x
  } else {
    lfx <- lf(x)
  }
  ly <- log(stunif()) + lfx
  # Step 2 ("Stepping out" procedure)
  L <- x - stunif() * w
  R <- L + w
  if ( ! is.finite(max_number_of_steps) ) {
    while ( ly < lf(L) ) {
      L <- L - w
    }
    while ( ly < lf(R) ) {
      R <- R + w
    }
  } else if ( max_number_of_steps > 0 ) {
    J <- floor( stunif() * max_number_of_steps )
    K <- ( max_number_of_steps - 1 ) - J
    while ( J > 0 && ly < lf(L) ) {
      L <- L - w
      J <- J - 1
    }
    while ( K > 0 && ly < lf(R) ) {
      R <- R + w
      K <- K - 1
    }
  }
  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    lfx1 <- lf(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, lfx=lfx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}

# Figure 5 of Neal (2003)
slice_sampler_shrinkage <- function(x, lf, L=0, R=1) {

  # Step 1
  if ( is.list(x) ) {
    lfx <- x$lfx
    x <- x$x
  } else {
    lfx <- lf(x)
  }
  ly <- log(stunif()) + lfx

  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    lfx1 <- lf(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, lfx=lfx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}

slice_sampler_latent <- function(x, s, lf, rate) {
  nEvaluations <- 1
  lfx <- lf(x)
  ly <- log(stunif()) + lfx
  half_s <- s/2
  l <- runif(1, x - half_s, x + half_s)
  # Eq. 7... a truncated exponential using the inverse CDF method.
  s <- -log(runif(1))/rate + 2*abs(l-x)
  half_s <- s/2
  L <- l - half_s
  R <- l + half_s
  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    nEvaluations <- nEvaluations + 1
    lfx1 <- lf(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, s=s, nEvaluations=nEvaluations)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}


# Algorithm of Li and Walker (2020)
slice_latent <- function(x, lf, lambda) {

  # Step 1
  if ( is.list(x) ) {
    lfx <- x$lfx
    s <- x$s
    x <- x$x
  } else {
    stop("Supplied x must be a list with previous value, previous latent s, and log-density eval at previous value.")
  }

  ly <- log(stunif()) + lfx

  # Step 2 sample latent variates
  ell <- runif(1, min = x - 0.5*s, max = x + 0.5*s)

  s_trunc <- 2.0*abs(ell - x)
  qtile_trunc <- 1.0 - exp(-lambda*s_trunc)
  s1 <- -log(1.0 - runif(1, min=qtile_trunc, max=1.0)) / lambda # draw from truncated exponential

  L <- ell - 0.5*s1
  R <- ell + 0.5*s1

  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    lfx1 <- lf(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, s=s1, lfx=lfx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}



## Algorithm 1 in Nishihara et al (2014)
slice_ellipse <- function(x, lf, mu, sig) {

  # Note that f is not the target density, but rather the target density divided by a normal prior density (with mean: mu, var: sig2),
  # i.e., the likelihood when the prior is normal

  # Step 1
  if ( is.list(x) ) {
    lfx <- x$lfx
    x <- x$x
  } else {
    stop("Supplied x must be a list with previous value and log-density eval at previous value.")
  }

  ly <- log(stunif()) + lfx

  # Step 2 (choose elipse)
  nu <- rnorm(1, mean=mu, sd=sig)
  theta <- runif(1, min=0.0, max=2.0*pi) # initial proposal

  L <- theta - 2.0*pi
  R <- theta

  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- (x - mu)*cos(theta) + (nu - mu)*sin(theta) + mu
    lfx1 <- lf(x1)
    if ( ly < lfx1 ) {
      result <- list(x=x1, lfx=lfx1)
      return(result)
    }
    if ( theta < 0.0 ) L <- theta else R <- theta
    theta <- runif(1, min=L, max=R)
  }
}

## Algorithm 2 in Nishihara et al (2014)
slice_ellipse_generalized <- function(x, lf, mu, sig, degf) {

  ## here, f is the target density

  a <- (degf + 1.0) / 2.0
  b <- 0.5*(degf + ((x$x - mu)/sig)^2)
  s <- 1.0 / rgamma(1, shape=a, rate=b) # rate of gamma <=> shape of inv-gamma

  lff <- function(xx) lf(xx) - ( dt((xx-mu)/sig, df=degf, log=TRUE) - log(sig) )

  slice_ellipse(x=x, lf=lff, mu=mu, sig=sqrt(s)*sig)
}
