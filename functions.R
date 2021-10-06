stunif <- function() runif(1, 0, 1)

# Figure 3 of Neal (2003)
slice_sampler_using_stepping_out <- function(x, f, w=1, max_number_of_steps=Inf) {
  # Step 1
  if ( is.list(x) ) {
    fx <- x$fx
    x <- x$x
  } else {
    fx <- f(x)
  }
  y <- stunif() * fx
  # Step 2 ("Stepping out" procedure)
  L <- x - stunif() * w
  R <- L + w
  if ( ! is.finite(max_number_of_steps) ) {
    while ( y < f(L) ) {
      L <- L - w
    }
    while ( y < f(R) ) {
      R <- R + w
    }
  } else if ( max_number_of_steps > 0 ) {
    J <- floor( stunif() * max_number_of_steps )
    K <- ( max_number_of_steps - 1 ) - J
    while ( J > 0 && y < f(L) ) {
      L <- L - w
      J <- J - 1
    }
    while ( K > 0 && y < f(R) ) {
      R <- R + w
      K <- K - 1
    }
  }
  # Step 3 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    fx1 <- f(x1)
    if ( y < fx1 ) {
      result <- list(x=x1, fx=fx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}

# Figure 5 of Neal (2003)
slice_sampler_shrinkage <- function(x, f, L=0, R=1) {
  
  # Step 1
  if ( is.list(x) ) {
    fx <- x$fx
    x <- x$x
  } else {
    fx <- f(x)
  }
  y <- stunif() * fx
  
  # Step 2 ("Shrinkage" procedure)
  repeat {
    x1 <- L + stunif() * ( R - L )
    fx1 <- f(x1)
    if ( y < fx1 ) {
      result <- list(x=x1, fx=fx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}



# Algorithm of Li and Walker (2020)
slice_latent <- function(x, f, lambda) {
  
  # Step 1
  if ( is.list(x) ) {
    fx <- x$fx
    s <- x$s
    x <- x$x
  } else {
    stop("Supplied x must be a list with previous value, previous latent s, and density eval at previous value.")
  }
  
  y <- stunif() * fx
  
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
    fx1 <- f(x1)
    if ( y < fx1 ) {
      result <- list(x=x1, s=s1, fx=fx1)
      return(result)
    }
    if ( x1 < x ) L <- x1 else R <- x1
  }
}


