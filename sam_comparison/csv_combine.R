library(tidyverse)
library(xtable)

fileNames <- list.files(path = '../data',
                        pattern = '.*.csv') %>% 
  paste0('../data/',.)

listOfDataFrames <- list(length = length(fileNames))

for(i in 1:length(fileNames)){
  listOfDataFrames[[i]] <- read_csv(fileNames[i]) %>% 
    cbind(curve = stringr::str_extract(fileNames[i],"curve[1-7]"))
}

df <- plyr::ldply(listOfDataFrames, data.frame) %>% 
  select(c('param_val','avgSampPSec','param','method','curve')) %>% 
  relocate('param', .before = 'param_val') %>% 
  relocate('avgSampPSec',.after = last_col())

curves <- stringr::str_extract(fileNames,"curve[1-7]") %>% unique()

for(i in 1:length(curves)){
tab <- df %>% filter(curve == curves[i]) %>% 
  rename('Parameter Value' = 'param_val',
         'Samples' = 'avgSampPSec',
         'Parameter' = 'param',
         'Method' = 'method',
         'Curve' = 'curve') %>% 
  xtable::xtable(caption = 'This table shows the parameter, parameter value, method of sampling, curve, and the average number of samples taken per second. Since these sampling methods lead to correlated samples the average number of samples per second was calculated after thinning the samples to remove auto correlation.',
                 digits = c(0,0,1,0,0,0)
                 )
xtable::print.xtable(tab, file = paste0(curves[i],'tbl.tex'))  
}

