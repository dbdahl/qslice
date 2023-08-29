# a function that takes the best set of parameters from each method
# author: Sam Johnson

combine_best_func <- function(list, metric = 'sampPsec') {
  gess_metrics <- list[[1]]
  latent_metrics <- list[[2]]
  rand_walk_metrics <- list[[3]]
  stepping_out_metrics <- list[[4]]
  transform_metrics <- list[[5]]
  # finding the best w
  temp <- stepping_out_metrics %>% select(w,eval(metric))
  temp <- do.call(rbind,lapply(split(temp, temp$w),FUN = colMeans)) %>% as.data.frame()
  if( metric == 'sampPsec' ) { 
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
  if( metric == 'sampPsec' ) {
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
  temp <- do.call(rbind,lapply(split(temp, temp$par),FUN = \(x) mean(x[,1]))) %>% as.data.frame()
  names(temp) <- 'sampPsec'
  if( metric == 'sampPsec' ) {
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
  if( metric == 'sampPsec' ) { 
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
  ##########
  # finding the best pseudo target
  temp <- transform_metrics %>% select(t,eval(metric)) |> mutate(t = unlist(t))
  temp <- do.call(rbind,lapply(split(temp, temp$t),FUN = \(list) mean(list$sampPsec))) %>% as.data.frame()
  temp$t <- row.names(temp)
  names(temp)[1] <- c('sampPsec')
  if( metric == 'sampPsec' ) { 
    bestPseudoTarget <- temp$t[which.max(temp[,metric])]
  } else {
    bestPseudoTarget <- temp$t[which.min(temp[,metric])]
  }
  ##########
  # creating a transform data frame
  tempTransform <- transform_metrics %>% 
    dplyr::filter(t == bestPseudoTarget) |> ### new
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
