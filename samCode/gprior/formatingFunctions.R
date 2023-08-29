## Function to help with formatting
# Author: Sam Johnson

attach(mtcars)
y <- scale(mpg)
X <- scale(cbind(cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb))
n <- length(y)
detach(mtcars)


# these functions are use to create a dataframe
singleLineFunc <- function(matrix, totTime) {
  vector <- c(matrix)
  effSamps <- coda::effectiveSize(vector)
  credInt <- round( quantile(vector, probs = c(0.025,0.975)), 2 )
  confInt <- round( mean(vector) + c(-1,0,1) * 1.96 * sd(vector)/sqrt(effSamps), 2 )
  SamplesPSec <- round(effSamps/totTime)
  c(credInt[1], lwrConf = confInt[1], est = confInt[2], uprConf = confInt[3], credInt[2],  SamplesPSec)
}


resultsFunc <- function(chainSamples) {
  # g parameter
  gSamples <- apply(chainSamples, 2, \(x) x$g) |> c()
  
  # psi parameter
  psiSamples <- apply(chainSamples, 2, \(x) x$psi) |> c()
  
  # beta parameters
  # betaSamples <- lapply(1:ncol(X), \(i) sapply(chainSamples, \(list) list$beta[,i]))
  
  # time
  time <- apply(chainSamples, 2, \(x) x$time) |> c()
  
  # list
  tempList <- list(gSamples) #, psiSamples)
  # tempList <- append(tempList, betaSamples)
  
  # formating the final Df
  tempList <- lapply(tempList, FUN = \(list) singleLineFunc(list, totTime = sum(time)))
  finalDf <- data.frame(do.call(rbind,tempList))
  names(finalDf) <- c('lwrCrd','lwrConf','Est','uprConf','uprCrd','SampPSec')
  row.names(finalDf) <- c('g')#,'psi',colnames(X))
  
  finalDf
}

# this creates traceplots
traceplot <- function(matrix, mainTitle = NULL) {
  n <- ncol(matrix)
  # finding eff sample size of samples
  effSize <- round(sum(apply(matrix, 2, effectiveSize)))
  # finding Rhat
  chainList <- mcmc.list()
  for (i in 1:n) chainList[[i]] <- mcmc(matrix[,i])
  Rhat <- round(gelman.diag(chainList)$psrf[1,1],3)
  # making trace plot
  color <- RColorBrewer::brewer.pal(n, 'Set2')
  plot(matrix[,1], type = 'l', col = color[i], main = paste0(mainTitle, ' ESS=',effSize, ' Rhat=',Rhat),
       ylab = '')
  for( i in 2:n ) {
    lines(matrix[,i], type = 'l', col = color[i])
  }
}