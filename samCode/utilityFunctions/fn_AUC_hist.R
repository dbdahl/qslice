auc <- function(x, y, u = NULL, nbins = 30) {

  if (is.null(u)) {
   
    if (is.function(y)) {
      y <- y(x)
    }
    
    stopifnot(length(x) == length(y))
    stopifnot(all(x > 0.0) && all(y >= 0.0) && all(x < 1.0))
    
  } else {

    bins = seq(0.0, 1.0, len = nbins + 1)
    x <- (bins[-(nbins+1)] + bins[-1]) / 2.0
    y <- tabulate( as.numeric(cut(u, breaks = bins)), nbins = nbins)
    
  }
  
  yn <- y / max(y)
  
  mean(yn)
}