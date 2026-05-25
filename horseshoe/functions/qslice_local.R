slice_quantile <- function (x, log_target, pseudo) 
{
  nEvaluations <- 0
  lf <- function(x) {
    nEvaluations <<- nEvaluations + 1
    log_target(x) - pseudo$ld(x)
  }
  
  lfx <- lf(x)
  
  if (is.nan(lfx) | is.infinite(lfx) | is.na(lfx)) {
    stop(paste0("input x = ", x, "\n",
                "log_target(x) = ", lfx, "\n"))
  }
  
  y <- log(runif(1, min = 0.0, max = 1.0)) + lfx
  
  L <- 0.0
  R <- 1.0
  
  repeat {
    
    u1 <- runif(1, L, R)
    x1 <- pseudo$q(u1)
    lfx1 <- lf(x1)
    
    if (is.nan(lfx1) | is.infinite(lfx1) | is.na(lfx1)) {
      
      if (x1 < pseudo$lb | x1 > pseudo$ub) {
        warning(paste0("Value proposed outside of boundary.\n",
                       "input x = ", x, "\n",
                       "proposed u1 = ", u1, "\n",
                       "proposed x1 = ", x1, "\n",
                       "log_target(x1) = ", lfx1, "\n",
                       "Proceeding with -Inf for lf(x1)"))
        lfx1 <- -Inf
      } else {
        warning(paste0("input x = ", x, "\n",
                       "proposed u1 = ", u1, "\n",
                       "proposed x1 = ", x1, "\n",
                       "log_target(x1) = ", lfx1, "\n"))
      }
      
    }
    
    if (y < lfx1) {
      return(list(x = x1, u = u1, nEvaluations = nEvaluations))
    }
    if (x1 < x) 
      L <- u1
    else R <- u1
  }
}