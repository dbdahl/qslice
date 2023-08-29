## water filling function
## author: Sam Johnson


# function to find how much water a curve can hold
water_area <- function(x, y, u = NULL, nbins = 30, plot = FALSE) {
  # browser()
  if(is.null(u)) {
    
    if (is.function(y)) {
      y <- y(x)
    } 
    
    stopifnot (length(y) == length(x))
    stopifnot(all(x > 0.0) && all(y >= 0.0) && all(x < 1.0))
    
  } else {
    
    bins <- seq(0.0, 1.0, len = nbins + 1)
    y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins=nbins)
    x <- (bins[-(nbins+1)] + bins[-1]) / 2.0
    
  }
  
  n <- length(y)
  
  # getting a normalized y
  yn <- y / max(y)
  # creating a vector to store which indecies are underwater
  underwater <- vector(length = n)
  # finding which values are underwater
  for( i in 2:(n-1) ) {
    if( max(yn[1:(i-1)]) > yn[i] && max(yn[(i+1):n]) > yn[i] ) {
      underwater[i] <- TRUE
    }
  }
  # checking if any points are underwater
  if( sum(underwater) == 0 ) {
    if( plot ) plot(x,yn,type = 'l')
    return(0)
  }
  # finding which indicies are under water
  underwaterIndex <- (1:n)[underwater]
  # finding the different bowls
  splits <- which(diff(underwaterIndex) != 1) + 1
  bowls <- unname(split(underwaterIndex, cumsum(seq_along(underwaterIndex) %in% splits)))
  # finding the areas of the bowls
  bowlArea <- sapply(1:length(bowls), FUN = \(i) {
    list <- bowls[[i]]
    waterLevel <- min(yn[c((min(list)-1),(max(list)+1))])
    area <- (waterLevel - yn[list]) |> sum()
    c(area = area, waterLevel = waterLevel)
  })
  # plotting the curve
  if( plot ) {
    plot(x,yn,bty = 'l', type = 'l', cex.axis = 1.75, cex.lab = 1.75)
    # points(x[underwater], yn[underwater], col = 'red')
    lapply(bowls, FUN = \(list) {
      for( i in 1:length(list) ) {
        index <- list[i]
        x0 <- x[index]
        y0 <- yn[index]
        y1 <- min(yn[c((min(list)-1),(max(list)+1))])
        segments(x0 = x0, x1 = x0, y0 = y0, y1 = y1, col = 'dodgerblue2')
      }
    })
  }
  sum(bowlArea[1,]) / n
}


# functions to test water fill function on
# fn <- \(x) dbeta(x,0.5,0.5); lf <- \(x) dbeta(x,0.5,0.5,log = TRUE) ;xx <- seq(from = .01, to = 1 - .01, length.out = 1000)
# fn <- \(x) 0.5 * dnorm(x) + 0.5 * dnorm(x,5,0.5); lf <- \(x) log(fn(x)) ;xx <- seq(from = -4, to = 9, length.out = 1000)
# fn <- \(x) 0.5 * dnorm(x) + 0.2 * dnorm(x,5,0.5) + 0.3 * dnorm(x,10,0.2); lf <- \(x) log(fn(x)); xx <- seq(from = -4, to = 12, length.out = 1000)
# fn <- \(x) 0.5 * dnorm(x) + 0.2 * dnorm(x,5,0.5) + 0.3 * dnorm(x,10,3); lf <- \(x) log(fn(x)); xx <- seq(from = -4, to = 20, length.out = 1000)
# fn <- \(x) 0.2 * dnorm(x) + 0.2 * dnorm(x,5,0.5) + 0.3 * dnorm(x,10,3) + 0.3 * dnorm(x,20,3.4); lf <- \(x) log(fn(x)); xx <- seq(from = -4, to = 30, length.out = 1000)
# fn <- \(x) sapply(x, \(x) ifelse(x < 5, 10-2*x, 2*x-10) * ifelse(x>=0 && x<=10,1,0)); lf <- \(x) log(fn(x)); xx <- seq(from = 0, to = 10, length.out = 1000)
# fn <- \(x) sapply(x, \(x) ifelse(x>=0 && x<=10, 20-2*x, 0)); lf <- \(x) log(fn(x)); xx <- seq(from = 0, to = 10, length.out = 1000)
# x <- seq(from = 0.01, to = 0.99, length = 20)
# y <- runif(length(x), 0.6, 1.0)
# water_area(x,y,plot = TRUE)
# 
# # creating the grid
# yy <- fn(xx)
# plot(xx,yy)
# 
# water_area(x = xx, y = yy, plot = TRUE)
# 
# # using samples to get the water function
# 
# nSamples <- 100000
# draws <- vector(length = nSamples)
# draws[1] <- mean(xx)
# for(i in 2:nSamples) {
#   draws[i] <- cucumber::slice_sampler_stepping_out(draws[i-1], target = lf, w = 100, max = Inf)$x
# } 
# 
# water_area(u = draws, plot = TRUE)
