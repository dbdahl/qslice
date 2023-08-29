# A function to format the data
# Author: Sam Johnson

format_df <- function(Df) {
  if( any(names(Df) %in% 'w') ) {
    Df <- Df %>%
      dplyr::select(nEval, sampPsec, w) %>% 
      rename('par' = 'w') %>% 
      mutate(par = as.character(par))
    Df$method <- 'SteppingOut'
  } else if ( any(names(Df) %in% 'rate') ) {
    Df <- Df %>% 
      dplyr::select(nEval, sampPsec, rate) %>% 
      rename('par' = 'rate') %>% 
      mutate(par = as.character(par))
    Df$method <- 'Latent'
  } else if ( any(names(Df) %in% 'mu') ) {
    Df <- Df %>% 
      dplyr::select(nEval, sampPsec, mu, sigma, df) %>% 
      mutate(par = paste0(mu, ',', sigma, ',', df)) %>% 
      dplyr::select(nEval, sampPsec, par) %>% 
      mutate(par = as.character(par))
    Df$method <- 'GESS'
  } else if ( any(names(Df) %in% 'c') ) {
    Df <- Df %>%
      dplyr::select(nEval, sampPsec, c) %>% 
      rename('par' = 'c')
    Df$method <- 'RandWalk'
  } else if ( any(names(Df) %in% 't') ) {
    Df <- Df %>% 
      dplyr::select(nEval, sampPsec, t) %>% 
      rename('par' = 't') %>% 
      mutate(par = unlist(par))
    Df$method <- 'Transform'
  } else {
    print('I do not know what went wrong')
  }
  Df
}