#' Water-Holding Area Above Curve (histogram)
#'
#' Function to find how much water a curve can hold above curve (between high points)
#' approximated by a histogram.
#' This area serves as a measure of convexity/multimodality.
#'
#' Author: Sam Johnson
#'
#' @param x Numeric vector of histogram locations. (Not used if \code{u} is supplied).
#' @param y Numeric vector of histogram heights OR function evaluating the curve
#' for a given value of \code{u}.(Not used if \code{u} is supplied).
#' @param u Numeric vector of samples supported on unit interval with which to
#' create histogram (use \code{u = NULL} if \code{x} and \code{y} are supplied).
#' @param nbins Number of histogram bins to use (defaults to 30).
#' @param plot Logical, whether to plot a visualization of the result.
#' @export
#' @examples
#' fn <- \(x) 0.5 * dbeta(x, 7, 2) + 0.5 * dbeta(x, 1, 6);
#' lf <- \(x) log(fn(x))
#' xx <- seq(from = 0.001, to = 0.999, length.out = 1000)
#' yy <- fn(xx)
#' water_area(x = xx, y = yy, plot = TRUE)
water_area <- function(x, y, u = NULL, nbins = 30, plot = FALSE) {
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
  # creating a vector to store which indices are underwater
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
