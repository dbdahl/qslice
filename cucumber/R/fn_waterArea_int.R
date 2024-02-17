#' Water-Holding Area Above Curve (function)
#'
#' Function to find how much water a curve can hold above curve (between high points).
#' This area serves as a measure of convexity/multimodality.
#'
#' Author: Sam Johnson
#'
#' @param f Univariate function (curve) for which to evaluate water-holding area.
#' @param interval Numeric vector of length two bounding the domain.
#' @param plot Logical, whether to plot a visualization of the result.
#' @param title String giving a title for the plot.
#' @param eps Positive numeric scalar for a small buffer on the interval. Defaults to \code{1.0e-3}.
#' @param tol_int Positive numeric scalar that passes to \code{abs.tol} in the call to \code{integrate()}.
#' Defaults to \code{1.0e-3}.
#' @returns List with:
#'
#' \code{totalArea}: Total area of a rectangle enclosing the curve (width is the
#' length of the interval, height is the max of \code{f} in the interval)
#'
#' \code{AUC}: Area under the curve (not restricted to unit square)
#'
#' \code{totalWaterArea}: Total water-holding area (not restricted to unit square)
#'
#' @export
#' @examples
#' f <- function(x) {0.33 * dnorm(x,0,0.5) + 0.33 * dnorm(x,2,1) + 0.33 * dnorm(x,7,2)}
#' interval <- c(-3,15)
#' water_area_int(f = f, interval = interval, plot = TRUE)
water_area_int <- function(f, interval = c(0,1), plot = FALSE, title = FALSE, eps = 1.0e-3, tol_int = 1.0e-3) {

  if (is.infinite(f(interval[1])) || is.nan(f(interval[1]))) {
    interval[1] <- interval[1] + eps
  }
  if (is.infinite(f(interval[2])) || is.nan(f(interval[2]))) {
    interval[2] <- interval[2] - eps
  }

  auc <- integrate(f, lower = interval[1], upper = interval[2], abs.tol = tol_int)$value
  df <- function(x) numDeriv::grad(f,x) # creating the derivative function
  extrema <- sort(rootSolve::uniroot.all(df, interval = interval + c(eps, -eps))) # finding all the extrema over the interval
  extrema <- c(interval[1], extrema, interval[2])

  nExtrema <- length(extrema) # finding the number of extrema

  f_at_extrema <- sapply(extrema, f)
  totalArea <- max(f_at_extrema) * abs(diff(interval))

  if (nExtrema > 2) {
    underwater <- sapply(2:(nExtrema - 1), \(i) {
      f_at_extrema[i] < max(f_at_extrema[i+1:(nExtrema-i)]) && f_at_extrema[i] < max(f_at_extrema[1:(i-1)]) # finding the extrema underwater
    })
  } else {
    underwater <- c(FALSE, FALSE)
  }

  if(sum(underwater) == 0) {
    if(plot) {
      main <- paste0('total area: ', round(totalArea,2), ' auc: ', round(auc,2), ' water area: ', 0)
      curve(f(x),interval[1],interval[2], ylim = c(0, max(f_at_extrema)), bty = 'l',
            main = ifelse(title,main,'')) # plotting the curve
      points(extrema, f(extrema), col = 'red') # plotting the extrema
    }
    return(list(totalArea = totalArea, AUC = auc, totalWaterArea = 0))
  }
  underwater <- c(FALSE, underwater, FALSE) # the edge points can't be underwater
  underwaterIndex <- (1:nExtrema)[underwater] # the index for each underwater extrema
  splits <- which(diff(underwaterIndex) != 1) + 1 # finding the bowls (non-contiguous underwater groups)
  bowls <- unname(split(underwaterIndex, cumsum(seq_along(underwaterIndex) %in% splits)))

  tst <- sapply(seq_along(bowls), \(i) {
      list <- bowls[[i]]
      minMax <- c(min(list)-1, max(list)+1) # finding the edges of the bowl
      extrema_now <- c(extrema[minMax]) # finding the x value at those edges
      waterLevel <- min(f(extrema_now)) # finding the value of that lower value

      tempfn <- \(x) f(x) - waterLevel # creating a temp function so I can find the range of the rectangle
      lwrUpr <- sort(rootSolve::uniroot.all(tempfn, interval = extrema_now)) # finding the upper and lower bound of the rectangle

      if(length(lwrUpr) < 2) { # in case uniroot missed the x at which it crosses (would happen if the interval is very small and therefore not important...)
        out <- list(waterArea = 0.0)
      } else {
        rectangleArea <- abs(diff(lwrUpr)) * waterLevel # finding the area of the rectangle
        areaUnderSection <-  integrate(f, lower = lwrUpr[1], upper = lwrUpr[2], abs.tol = tol_int) # finding the area under the curve
        waterArea <- rectangleArea - areaUnderSection$value # the area the water is filling
        out <- list(waterLevel = waterLevel, lwr = lwrUpr[1], upr = lwrUpr[2], rectangleArea = rectangleArea, waterArea = waterArea, areaUnderSection = areaUnderSection$value)
      }
    }, simplify = FALSE
  )

  totalWaterArea <- sum(sapply(seq_along(tst), \(i) {tst[[i]]$waterArea}))

  if(plot) {
    main <- paste0('total area: ', round(totalArea,2), ' auc: ', round(auc,2), ' water area: ', round(totalWaterArea,2))
    curve(f(x),interval[1],interval[2], ylim = c(0, max(f_at_extrema)), bty = 'l',
          main = ifelse(title,main,'')) # plotting the curve
    points(extrema, f(extrema), col = 'red') # plotting the extrema
    lapply(tst, \(list) { # creating the water lines
      segments(x0 = list$upr, x1 = list$lwr, y0 = list$waterLevel, y1 = list$waterLevel, col = 'blue')
      text(x = mean(c(list$upr, list$lwr)), y = list$waterLevel*0.960, paste0('water area: ', round(list$waterArea,3)),
           cex = 0.75, font = 3, col = 'blue')
    })
  }
  list(totalArea = totalArea, AUC = auc, totalWaterArea = totalWaterArea)
}
