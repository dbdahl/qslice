# this is a summary of the results from the samplers with a normal target
# author: Sam Johnson

# sourcing the files the different samplers

library(doParallel)
library(parallel)
library(tidyverse)

auc_diagnostic = function(samples_u, nbins = 30) {
  
  (bins = seq(0.0, 1.0, len=nbins+1))
  (tab = tabulate( as.numeric(cut(samples_u, breaks = bins)), nbins=nbins))
  
  tab[c(1,nbins)] = 1.2 * tab[c(1,nbins)] # penalty for smiling...
  
  (tab_norm = tab / max(tab) / nbins)
  
  (auc = sum(tab_norm))
  auc
}

files <- Sys.glob('data/*.rds')
output <- sapply(files, FUN = readRDS)
targetTitle <- 'Normal Target'

transform_metrics <- output$`data/transform.rds`
stepping_out_metrics <- output$`data/stepping_out.rds`
latent_metrics <- output$`data/latent.rds`
gess_metrics <- output$`data/gess.rds`
rand_walk_metrics <- output$`data/randWalk.rds`

# summary plot
## transform
pdf(file = 'images/transform_nEval.pdf')
transform_metrics %>% 
  ggplot(aes(y = unlist(t), x = nEval)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'psuedo target'
  )
dev.off()

pdf(file = 'images/transform_SampPSec.pdf')
transform_metrics %>% 
  ggplot(aes(y = unlist(t), x = SampPSec)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'psuedo target'
  )
dev.off()

pdf(file = 'images/transform_SampPSecAndnEval.pdf')
transform_metrics %>% 
  ggplot(aes(x = nEval, y = SampPSec, color = unlist(t))) + 
  geom_point(alpha = 0.5) + 
  geom_point(data = transform_metrics %>% 
               group_by(unlist(t)) %>% 
               summarise(nEval = mean(nEval), SampPSec = mean(SampPSec)) %>% 
               rename('t' = 'unlist(t)'),
             aes(x = nEval, y = SampPSec, color = t), size = 5, shape = 9)
dev.off()                

## stepping out
pdf(file = 'images/steppingOut_nEval.pdf')
stepping_out_metrics %>% 
  ggplot(aes(y = factor(w), x = nEval)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'w'
  )
dev.off()

pdf(file = 'images/steppingOut_SampPSec.pdf')
stepping_out_metrics %>% 
  ggplot(aes(y = factor(w), x = SampPSec)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'w'
  )
dev.off()

pdf(file = 'images/steppingOut_SampPSecAndnEval.pdf')
stepping_out_metrics %>% 
  dplyr::filter(w > 3) %>% 
  ggplot(aes(x = nEval, y = SampPSec, color = factor(w))) +
  geom_point(alpha = 0.5) +
  geom_point(data = stepping_out_metrics %>% 
               dplyr::filter(w > 3) %>% 
               group_by(w) %>% 
               summarise(SampPSec = mean(SampPSec), nEval = mean(nEval)),
             aes(x = nEval, y = SampPSec, color = factor(w)), size = 5, shape = 9)
dev.off()

## latent
pdf(file = 'images/latent_nEval.pdf')
latent_metrics %>% 
  ggplot(aes(y = factor(rate), x = nEval)) +
  geom_boxplot() + 
  labs(
    title = targetTitle,
    y = 'rate'
  )
dev.off()

pdf(file = 'images/latent_SampPSec.pdf')
latent_metrics %>% 
  ggplot(aes(y = factor(rate), x = SampPSec)) +
  geom_boxplot() + 
  labs(
    title = targetTitle,
    y = 'rate'
  )
dev.off()

pdf(file = 'images/latent_SampPSecAndnEval.pdf')
latent_metrics %>% 
  ggplot(aes(x = nEval, y = SampPSec, color = factor(rate))) +
  geom_point(alpha = 0.5) +
  geom_point(data = latent_metrics %>% 
               group_by(rate) %>% 
               summarise(SampPSec = mean(SampPSec), nEval = mean(nEval)),
             aes(x = nEval, y = SampPSec, color = factor(rate)), size = 5, shape = 9)
dev.off()

## gess
pdf(file = 'images/gess_nEval.pdf')
gess_metrics %>% 
  ggplot(aes(y = paste0(mu,',',sigma,',',df), x = nEval)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'mu,sigma,df'
  )
dev.off()

pdf(file = 'images/gess_SampPSec.pdf')
gess_metrics %>% 
  ggplot(aes(y = paste0(mu,',',sigma,',',df), x = SampPSec)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'mu,sigma,df'
  )
dev.off()

pdf(file = 'images/gess_SampPSecAndnEval.pdf')
gess_metrics %>%
  mutate(par = paste0(mu,',',sigma,',',df)) %>% 
  ggplot(aes(x = nEval, y = SampPSec, color = factor(par))) +
  geom_point(alpha = 0.5) +
  geom_point(data = gess_metrics %>% 
               mutate(par = paste0(mu,',',sigma,',',df)) %>% 
               group_by(par) %>% 
               summarise(SampPSec = mean(SampPSec), nEval = mean(nEval)),
             aes(x = nEval, y = SampPSec, color = factor(par)), size = 5, shape = 9)
dev.off()

## rand walk
pdf(file = 'images/randWalk_AcceptanceRate.pdf')
rand_walk_metrics %>%
  ggplot(aes(y = factor(c), x = acceptanceRate)) +
  geom_boxplot() +
  labs(
    y = 'c'
  )
dev.off()

pdf(file = 'images/randWalk_SampPSec.pdf')
rand_walk_metrics %>% 
  ggplot(aes(y = factor(c), x = SampPSec)) +
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'c'
  )
dev.off()

# combining them all
combineFunc <- function(stepping_out_metrics, latent_metrics, gess_metrics, rand_walk_metrics, transform_metrics, metric) {
  # finding the best w
  temp <- stepping_out_metrics %>% select(w,eval(metric))
  temp <- do.call(rbind,lapply(split(temp, temp$w),FUN = colMeans)) %>% as.data.frame()
  if( metric == 'SampPSec' ) { 
    bestW <- temp$w[which.max(temp[,metric])]
  } else {
    bestW <- temp$w[which.min(temp[,metric])]
  }
  # creating a temporary stepping out data frame
  tempSteppingOut <- stepping_out_metrics %>% 
    dplyr::filter(w == bestW) %>% select(eval(metric),w) %>% 
    dplyr::rename('par' = 'w') %>% 
    mutate(par = paste0('w=',par))
  tempSteppingOut$method <- 'SteppingOut'
  # finding the best rate
  temp <- latent_metrics %>% select(rate,eval(metric))
  temp <- do.call(rbind,lapply(split(temp, temp$rate),FUN = colMeans)) %>% as.data.frame()
  if( metric == 'SampPSec' ) {
    bestRate <- temp$rate[which.max(temp[,metric])]
  } else {
    bestRate <- temp$rate[which.min(temp[,metric])]
  }
  # creating a temporary latent data frame
  tempLatent <- latent_metrics %>% 
    dplyr::filter(rate == bestRate) %>% select(eval(metric), rate) %>% 
    dplyr::rename('par' = 'rate') %>% 
    mutate(par = paste0('rate=',par))
  tempLatent$method <- 'Latent'
  # finding the best gess metrics
  temp <- gess_metrics %>% select(mu,sigma,df,eval(metric)) %>% mutate(par = paste0(mu,sigma,df)) %>% select(-c(mu,sigma,df))
  temp <- do.call(rbind,lapply(split(temp, temp$par),FUN = \(list) colMeans(list[,1]))) %>% as.data.frame()
  if( metric == 'SampPSec' ) {
    bestGESS <- rownames(temp)[which.max(temp[,metric])]
  } else {
    bestGESS <- rownames(temp)[which.min(temp[,metric])]
  }
  #  creating a temporary gess data frame
  tempGess <- gess_metrics %>% 
    mutate(par = paste0(mu,sigma,df)) %>% 
    dplyr::filter(par == bestGESS) %>% 
    mutate(par = paste0('mu=',mu,',sigma=',sigma, ',df=',df)) %>% 
    dplyr::select(eval(metric), par)
  tempGess$method <- 'GESS'
  # finding a the best c
  temp <- rand_walk_metrics %>% select(c,eval(metric))
  temp <- do.call(rbind,lapply(split(temp, temp$c),FUN = colMeans)) %>% as.data.frame()
  if( metric == 'SampPSec' ) { 
    bestC <- temp$c[which.max(temp[,metric])]
  } else {
    bestC <- temp$c[which.min(temp[,metric])]
  }
  # creating a temporary rand walk data frame
  tempRandWalk <- rand_walk_metrics %>% 
    dplyr::filter(c == bestC) %>% 
    select(eval(metric), c) %>% 
    mutate(c = paste0('c=',c)) %>% 
    dplyr::rename('par' = 'c')
  tempRandWalk$method <- 'RandWalk'
  # creating a transform data frame
  tempTransform <- transform_metrics %>% 
    select(eval(metric), t) %>% 
    mutate(t = unlist(t)) %>% 
    mutate(
      t = stringr::str_extract(t, pattern = 'l.*') %>%
        gsub('[)]','',.)
    ) %>% 
    dplyr::rename('par' = 't')
  tempTransform$method <- 'Transform'
  # combining all the data
  allSamplers <- do.call(rbind, list(tempTransform,tempSteppingOut, tempLatent, tempGess, tempRandWalk))
  
  allSamplers
}

totalSampPSec <- combineFunc(stepping_out_metrics = stepping_out_metrics, 
            latent_metrics = latent_metrics,
            gess_metrics = gess_metrics,
            rand_walk_metrics = rand_walk_metrics,
            transform_metrics = transform_metrics,
            metric = 'SampPSec')

pdf(file = 'images/totalSampPSec.pdf')
totalSampPSec %>% 
  ggplot(aes(y = par, x = SampPSec, fill = method)) + 
  geom_boxplot() +
  labs(
    title = targetTitle
  )
dev.off()

totalNEval <- combineFunc(stepping_out_metrics = stepping_out_metrics, 
            latent_metrics = latent_metrics,
            gess_metrics = gess_metrics,
            rand_walk_metrics = rand_walk_metrics,
            transform_metrics = transform_metrics,
            metric = 'nEval')

pdf(file = 'images/totalNEval.pdf')
totalNEval %>% 
  ggplot(aes(y = par, x = nEval/50000, fill = method)) + 
  geom_boxplot() +
  labs(
    title = targetTitle
  )
dev.off()

formatDf <- function(Df) {
  if( any(names(Df) %in% 'w') ) {
    Df <- Df %>%
      dplyr::select(nEval, SampPSec, w) %>% 
      rename('par' = 'w') %>% 
      mutate(par = as.character(par))
    Df$method <- 'SteppingOut'
  } else if ( any(names(Df) %in% 'rate') ) {
    Df <- Df %>% 
      dplyr::select(nEval, SampPSec, rate) %>% 
      rename('par' = 'rate') %>% 
      mutate(par = as.character(par))
    Df$method <- 'Latent'
  } else if ( any(names(Df) %in% 'mu') ) {
    Df <- Df %>% 
      dplyr::select(nEval, SampPSec, mu, sigma, df) %>% 
      mutate(par = paste0(mu, ',', sigma, ',', df)) %>% 
      dplyr::select(nEval, SampPSec, par) %>% 
      mutate(par = as.character(par))
    Df$method <- 'GESS'
  } else if ( any(names(Df) %in% 'c') ) {
    Df <- Df %>%
      dplyr::select(nEval, SampPSec, c) %>% 
      rename('par' = 'c')
    Df$method <- 'RandWalk'
  } else if ( any(names(Df) %in% 't') ) {
    Df <- Df %>% 
      dplyr::select(nEval, SampPSec, t) %>% 
      rename('par' = 't') %>% 
      mutate(par = unlist(par))
    Df$method <- 'Transform'
  } else {
    print('I do not know what went wrong')
  }
  Df
}

total_metrics <- lapply(list(stepping_out_metrics, latent_metrics, gess_metrics, rand_walk_metrics, transform_metrics),
                        formatDf)
total_metrics <- do.call(rbind, total_metrics)

pdf(file = 'images/total_SampPSecAndnEval.pdf')
total_metrics %>%
  ggplot(aes(x = nEval, y = SampPSec, color = method)) +
  geom_point(alpha = 0.15) +
  geom_point(data = total_metrics %>% 
               group_by(par, method) %>% 
               summarise(SampPSec = mean(SampPSec), nEval = mean(nEval)),
             aes(x = nEval, y = SampPSec, color = method), size = 5, shape = 9)
dev.off()

## ploting what the transformations look like

avgSampPSec <- totalSampPSec %>% 
  group_by(par) %>% 
  summarise(SampPSec = mean(SampPSec)) 

compDf <- avgSampPSec %>% 
  rowwise() %>% 
  mutate(bStep = ifelse(SampPSec >= avgSampPSec[grepl('.*w.*',avgSampPSec$par),"SampPSec"], 1, 0),
         bGess = ifelse(SampPSec >= avgSampPSec[grepl('mu=.*',avgSampPSec$par),"SampPSec"], 1, 0),
         bLatent = ifelse(SampPSec >= avgSampPSec[grepl('rate=.*',avgSampPSec$par),"SampPSec"], 1, 0),
         bRandWalk = ifelse(SampPSec >= avgSampPSec[grepl('^c=[1-9]*',avgSampPSec$par),"SampPSec"], 1, 0),
         bTotal = sum(bStep, bGess, bLatent, bRandWalk)
           ) %>% 
  slice(grep("^loc.*", par))


uniqueIndex <- !duplicated(transform_metrics$t)
pseu <- vector('list', length = sum(uniqueIndex))
uniqueT <- transform_metrics$t[uniqueIndex]
uniquePDF <- transform_metrics$log_pdf[uniqueIndex]
uniqueInvCDF <- transform_metrics$inv_cdf[uniqueIndex]
colFunc <- colorRampPalette(c('red','blue'))
colors <- RColorBrewer::brewer.pal(4, 'Dark2')

temp <- unlist(uniqueT)
temp <- stringr::str_extract(temp, pattern = 'l.*') %>%
  gsub('[)]','',.)


compDf <- compDf[match(temp, compDf$par),]

for( i in seq_along(pseu) ) {
  pseu[[i]] <- list(d = function(x) exp(uniquePDF[[i]](x)),
                    q = function(u) uniqueInvCDF[[i]](u),
                    t = uniqueT[[i]],
                    color = colors[compDf$bTotal[i]],
                    lwd = compDf$bTotal[i],
                    lty = compDf$bTotal[i])
}

truth <- list( d = function(x) dnorm(x, 0, 1),
               q = function(u) qnorm(u, 0, 1),
               t = 'norm(0,1)')

### plot
(n_pseu = length(pseu))

xx <- seq(from = -4, to = 4, length.out = 1000)
uu <- seq(from = 0, to = 1, length.out = 1000)

pdf(file = 'images/plotOfTransformation.pdf')
m = matrix(c(1,2,3), nrow=1)
layout(mat=m, widths=c(1, 1, 0.7))
par(mar=c(4,1,1,1))

plot(xx, truth$d(xx), type="l", lwd=5,
     ylim=c(0,.6),
     xlab=expression(theta), ylab="", axes=FALSE)
axis(side=1)
for(i in 1:n_pseu) {
  lines(xx, pseu[[i]]$d(xx), col = pseu[[i]]$color, lwd = pseu[[i]]$lwd, lty = pseu[[i]]$lty)#col=i+1, lwd=2, lty=i+1)	
}

plot(uu, rep(1.0, length(uu)), type="l", col=1, lwd=5, 
     ylim=c(0,10),
     xlab=expression(psi), ylab="", axes=FALSE, bty="L")
axis(side=1)

for(i in 1:n_pseu) {
  lines(uu, truth$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col = pseu[[i]]$color, lwd = pseu[[i]]$lwd, lty = pseu[[i]]$lty)#col=i+1, lwd=2, lty=i+1)
}

plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("left", bty="n", inset=c(0.0, 0.0),
       legend = 1:4,
       col = colors, lwd = 1:4, cex = 0.9, lty = 1:4)
       # legend=c(truth$t, sapply(pseu, function(x) x$t)),
       # col=1:(n_pseu+1), lwd=2, lty=1:(n_pseu+1), cex = 0.7)
dev.off()



# plotting the us

transformUdraws <- transform_metrics |> 
  dplyr::select(t, udraws) |>
  mutate(t = unlist(t)) |> 
  group_by(t) |> 
  reframe(udraws = list(unlist(list(udraws))))

pdf(file = 'images/histOfUDraws.pdf')
par(mfrow = c(ceiling(nrow(transformUdraws)/3),3))

apply(transformUdraws, 1, FUN = \(row) {
  psuedoTarget <- row[1]
  u <- unlist(row[2])
  aucDiagnostic <- round(auc_diagnostic(u),3)
  hist(u, main = psuedoTarget$t, sub = paste0('auc:',aucDiagnostic))
})
dev.off()