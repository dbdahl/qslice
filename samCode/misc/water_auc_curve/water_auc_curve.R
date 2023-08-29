## testing how water and auc effect the slice sampler
## author: sam johnson

# loading in neccesary libraries
library(tidyverse)
library(plotly)

# loading necessary scripts to calculate auc and water

source("../utility_functions/fn_AUC_hist.R")
source("../utility_functions/fn_skew_hist.R")
source("../utility_functions/fn_waterArea_hist.R")
source("../utility_functions/fn_waterArea_int.R")
source("../utility_functions/fn_expectedSliceWidth.R")
source("../utility_functions/fn_shrinkslice_utility.R")


# slice sampler that only uses shrinkage
slice_sampler_shrinkage <- function (x, target, interval = c(0,1), log = TRUE) 
{
  nEvaluations <- 0
  f <- function(x) {
    nEvaluations <<- nEvaluations + 1
    target(x)
  }
  fx <- f(x)
  y <- if (isTRUE(log)) {
    log(runif(1)) + fx
  }
  else {
    runif(1) * fx
  }
  L <- interval[1]
  R <- interval[2]
  
  repeat {
    x1 <- L + runif(1) * (R - L)
    if (y < f(x1)) 
      return(list(x = x1, nEvaluations = nEvaluations))
    if (x1 < x) 
      L <- x1
    else R <- x1
  }
}

shrinkage_eval <- function(samples, int.x, lf, interval) {
  draws <- vector(length = samples); nEval <- vector(length = samples)
  draws[1] <- int.x
  time <- system.time({
    for(i in 2:samples) {
      temp <- slice_sampler_shrinkage(draws[i-1], target = lf, interval = interval, log = TRUE)
      draws[i] <- temp$x
      nEval[i] <- temp$nEvaluations
    }
  })
  util <- utility_shrinkslice(h = \(x) exp(lf(x)), type = 'function')
  effSamps <- coda::effectiveSize(draws)
  list(sampsPSec = effSamps/time[1], samples = samples, effSamps = coda::effectiveSize(draws),
       time = time, draws = draws, nEval = nEval, util = util)
}

# function to test
x = seq(0.01, 0.99, length=100)
dev.off()
par(mfrow=c(1,2))
y = \(x) dbeta(x, 0.5, 0.5); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) dbeta(x, 0.5, 0.5, log = TRUE)
y = \(x) dbeta(x, 1.0, 0.8); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) dbeta(x, 1.0, 0.8, log = TRUE) 
y = \(x) dbeta(x, 1.0, 1.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) dbeta(x, 1.0, 1.0, log = TRUE)
y = \(x) dbeta(x, 2.0, 2.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) dbeta(x, 2.0, 2.0, log = TRUE)
y = \(x) dbeta(x, 1.0, 2.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) dbeta(x, 1.0, 2.0, log = TRUE)
y = \(x) dbeta(x, 20.0, 20.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) dbeta(x, 20.0, 20.0, log = TRUE)
y = \(x) dbeta(x, 3.0, 8.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE);  lf <- \(x) dbeta(x, 3.0, 8.0, log = TRUE)
y = \(x) 0.9*dbeta(x, 3.0, 8.0) + 0.1*dbeta(x, 15, 1.5); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) log(0.9*dbeta(x, 3.0, 8.0) + 0.1*dbeta(x, 15, 1.5))
y = \(x) 0.5*(x >= 0.5)*dbeta(x, 1.0, 2.0) + 0.5*(x < 0.5)*dbeta(x, 2.0, 1.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) log(0.5*(x >= 0.5)*dbeta(x, 1.0, 2.0) + 0.5*(x < 0.5)*dbeta(x, 2.0, 1.0))
y = \(x) 0.5*(x < 0.5)*(1.0 - 2.0*x) + 0.5*(x >= 0.5)*2.0*(x - 0.5); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) log(0.5*(x < 0.5)*(1.0 - 2.0*x) + 0.5*(x >= 0.5)*2.0*(x - 0.5))
y = \(x) 0.4*dbeta(x, 30.0, 10.0) + 0.6*dbeta(x, 10.0, 30.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) log(0.4*dbeta(x, 30.0, 10.0) + 0.6*dbeta(x, 10.0, 30.0))
y = \(x) 0.4*dbeta(x, 30.0, 3.0) + 0.6*dbeta(x, 3.0, 30.0) + 0.2; utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) log(0.4*dbeta(x, 30.0, 3.0) + 0.6*dbeta(x, 3.0, 30.0) + 0.2)
y = \(x) 2.0*dbeta(x, 30.0, 10.0) + 5.0*dbeta(x, 3, 3) + 1.7*dbeta(x, 10.0, 30.0); utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE); lf <- \(x) log(2.0*dbeta(x, 30.0, 10.0) + 5.0*dbeta(x, 3, 3) + 1.7*dbeta(x, 10.0, 30.0))
y = \(x) 2.0*dbeta(x,9,36) + 5.0*dbeta(x,30,30) + 2.7*dbeta(x,50,12);  utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE); utility_shrinkslice(h = y, type = "function", plot = TRUE)
y = runif(length(x), in = 0.6, max = 1.0); utility_shrinkslice(x = x, y = y, type = "grid", plot = TRUE)
dev.off()

pdf(file = 'exampleOfWater.pdf')
# utility_shrinkslice(h = y, type = "grid", plot = TRUE) # if using the function
utility_shrinkslice(x = x, y = y(x), type = "grid", plot = TRUE)
dev.off()

functionsToTest <- list(
  \(x) dbeta(x, 0.5, 0.5, log = TRUE),
  \(x) dbeta(x, 1.0, 0.8, log = TRUE),
  \(x) dbeta(x, 1.0, 1.0, log = TRUE),
  \(x) dbeta(x, 2.0, 2.0, log = TRUE),
  \(x) dbeta(x, 1.0, 2.0, log = TRUE),
  \(x) dbeta(x, 20.0, 20.0, log = TRUE),
  \(x) dbeta(x, 3.0, 8.0, log = TRUE),
  \(x) log(0.9*dbeta(x, 3.0, 8.0) + 0.1*dbeta(x, 15, 1.5)),
  \(x) log(0.5*(x >= 0.5)*dbeta(x, 1.0, 2.0) + 0.5*(x < 0.5)*dbeta(x, 2.0, 1.0)),
  \(x) log(0.5*(x < 0.5)*(1.0 - 2.0*x) + 0.5*(x >= 0.5)*2.0*(x - 0.5)),
  \(x) log(0.4*dbeta(x, 30.0, 10.0) + 0.6*dbeta(x, 10.0, 30.0)),
  \(x) log(0.4*dbeta(x, 30.0, 3.0) + 0.6*dbeta(x, 3.0, 30.0) + 0.2),
  \(x) log(2.0*dbeta(x, 30.0, 10.0) + 5.0*dbeta(x, 3, 3) + 1.7*dbeta(x, 10.0, 30.0))
)

# test
samplesFromCurves <- sapply(functionsToTest, FUN = \(x) {
  functionName <- deparse(x) |>
    stringr::str_remove_all('log|, log = TRUE') |>
    stringr::str_remove_all('x,') |>
    stringr::str_replace_all(' ','') #|>
    # stringr::str_replace_all('([[(]])\\1+', '\\1')
    
  # cant do as.character(x) found way to convert a function to a string
  temp <- replicate(10, shrinkage_eval(5e4, int.x = 0.5, lf = x, interval = c(0,1)))
  temp <- apply(temp, 2, \(x) {
    data.frame(sampsPSec = x$sampsPSec, util = x$util[1], auc = x$util[2], water = x$util[3])
  })
  temp <- do.call(rbind, temp)
  temp$curve <- functionName[2]
  temp
}, simplify = FALSE)

resultsDF <- do.call(rbind, samplesFromCurves)
 
save(resultsDF, file = 'resultsDF.rda')

load(file = 'resultsDF.rda')

resultsDF |> 
  rowwise() |> 
  mutate(auc = auc + rnorm(1, 0, 0.01),
         water = water + rnorm(1, 0, 0.01)) |> 
  ggplot(aes(x = auc, y = water, color = sampsPSec)) +
  geom_point()


fig <- plot_ly(resultsDF, x = ~auc, y = ~water, z = ~sampsPSec, color = ~curve, colors = RColorBrewer::brewer.pal(13,'Paired'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'AUC'),
                                   yaxis = list(title = 'Water'),
                                   zaxis = list(title = 'Samples Per Second')))

fig
