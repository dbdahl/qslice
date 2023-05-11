## This scripts creates all the graphics ##
## author: Sam Johnson ##

# sourcing needed files
source('sampler_eda.R')
# choosing the folder to save the images to
file_path <- '../images_slice_sampler_comp/sampler_eda/'
save.images <- FALSE 
heightVal <-  1.445
widthVal <- 2.83

curveNames <- c('norm(20,5)', 'gamma(2.5,1)',
                '.2 norm(0,0.5) + .8 norm(6,2)',
                '.2 norm(0,0.5) + .8 norm(20,1)',
                't(5)','truncated t(5)','beta(0.2,0.8)')

# this plots each curve
if(save.images) {pdf(file = paste0(file_path, 'curves.pdf'))}
par(mfrow = c(3,3))
for (i in 1:(length(curves) - 1)) {
  curve(fexp(f = curves[[i]], x = x), 
        xlim = range[[i]]$xrange,
        ylim = range[[i]]$yrange,
        ylab = "", xlab = '', main = curveNames[i])
}
if(save.images) {dev.off()}


# reading in all files from data folder
files <- list.files(path = '../data/', pattern = "*_metrics", full.names = TRUE)
# saving the files to a dataframe
curve_comp_data <- lapply(files, extract_metrics) %>% plyr::ldply(data.frame)

names(curveNames) <- c('7','3','1','2','5','4','6')

# Overall Performance
if(save.images) pdf(file = paste0(file_path, 'curvesPerformance.pdf'))
curve_comp_data %>% 
  dplyr::filter(curve != 8) %>% 
  dplyr::filter(ksTest > 0.01) %>% 
  mutate(method = case_when(
    method == 'transform' ~ 'Transform',
    method == 'stepping_out' ~ 'Neal Stepping Out',
    method == 'rand_walk' ~ 'Random Walk',
    method == 'latent' ~ 'Latent',
    method == 'gess' ~ 'GESS',
    TRUE ~ method
  ),
  curve = factor(curve, levels = c('7','3','1','2','5','4','6'))) %>% 
  ggplot(aes(y = method, x = SampPSec, fill = method)) +
  facet_wrap(~curve, labeller = labeller(curve = curveNames)) + 
  geom_boxplot() +
  scale_x_continuous(labels = scales::number_format(scale = .001, suffix = 'K')) +
  labs(
    x = 'Effective Samples Per CPU Second',
    y = 'Method'
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 8)
  )
if(save.images) dev.off()
rm(curve_comp_data)


# i for each curve

for(i in 1:7) {
  print(i)
  # this plots the curve by itself
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'.pdf'))
  curve(fexp(f = curves[[i]], x = x), xlim = range[[i]]$xrange, ylim = range[[i]]$yrange, main = curveNames[i], xlab = '', ylab = '')
  if(save.images) dev.off()
 
  ksFilterInd <- TRUE

  # this extracts all the draws
  draws_list <- lapply(files[str_detect(files, paste0('.*',i,'.*_metrics'))], extract_draws)
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'Results.pdf'))
  # plots the density of all the curves against the true curve
  curve_slice_sampler(draws_list,curve = i, ksFilter = FALSE)
  if(save.images) dev.off()
  # removes the excess data
  rm(draws_list)
  gc()
  
  
  
  # extracts the data about each curve except the draws
  method_comp_list <- lapply(files[str_detect(files, paste0('.*',i,'.*_metrics'))], extract_curve_metrics)
  # this function plots all the boxplots
  curve <- boxplot_slice_sampler(method_comp_list, curve = i, ksFilter = ksFilterInd)
  # removes extra data
  rm(method_comp_list)
  gc()
  
  # a plot comparing each slice sampler
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'SamplerComp.pdf'))
  print(curve$sampler_comp)
  if(save.images) dev.off()
  
  
  ### Stepping Out
  method <- 'stepping_out'
  
  # boxplot of x values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'SteppingoutXBoxplot.pdf'))
  print(curve$stepping_out_x)
  if(save.images) dev.off()
  # boxplot of w values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'SteppingoutWBoxplot.pdf'))
  print(curve$stepping_out_w)
  if(save.images) dev.off()
  # boxplot of parameters
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'SteppingoutParamBoxplot.pdf'))
  print(curve$stepping_out_param)
  if(save.images) dev.off()
  # boxplot of thining
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'SteppingoutThinBoxplot.pdf'))
  print(curve$stepping_out_samplesThin)
  if(save.images) dev.off()
  
  
  # extracting the draws for the stepping out procedure
  draws_list <- lapply(files[str_detect(files, paste0('.*',i,'.*', method,'_metrics'))], extract_draws)
  # plot the density of the draws against true value
  if(save.images) pdf(file = paste0(file_path, 'curve',i, method, 'Density.pdf'))
  facet_plot(draws_list)
  if(save.images) dev.off()
  # removing excess data
  rm(draws_list)
  gc()
  
  
  
  ### Latent
  method <- 'latent'
  
  # boxplot of x values 
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'LatentXBoxplot.pdf'))
  print(curve$latent_x)
  if(save.images) dev.off()
  # boxplot of s values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'LatentSBoxplot.pdf'))
  print(curve$latent_s)
  if(save.images) dev.off()
  # boxplot of rate values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'LatentRateBoxplot.pdf'))
  print(curve$latent_rate)
  if(save.images) dev.off()
  # boxplot of parameters
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'LatentParamBoxplot.pdf'))
  print(curve$latent_param)
  if(save.images) dev.off()
  # boxplot of thinning
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'LatentThinBoxplot.pdf'))
  print(curve$latent_samplesThin)
  if(save.images) dev.off()
  
  
  # extracting the draws for the latent
  draws_list <- lapply(files[str_detect(files, paste0('.*',i,'.*', method,'_metrics'))], extract_draws)
  # plots the density of the draws
  if(save.images) pdf(file = paste0(file_path, 'curve',i, method,'Density.pdf'))
  facet_plot(draws_list)
  if(save.images) dev.off()
  # removing excess data
  rm(draws_list)
  gc()
  
  
  
  ### GESS
  method <- 'gess'
  
  # boxplot of x values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'GessXBoxplot.pdf'))
  print(curve$gess_x)
  if(save.images) dev.off()
  # boxplot of mu values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'GessMuBoxplot.pdf'))
  print(curve$gess_mu)
  if(save.images) dev.off()
  # boxplot of sigma values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'GessSigmaBoxplot.pdf'))
  print(curve$gess_sigma)
  if(save.images) dev.off()
  # boxplot of df values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'GessDfBoxplot.pdf'))
  print(curve$gess_df)
  if(save.images) dev.off()
  # boxplot of parameter values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'GessParamBoxplot.pdf'))
  print(curve$gess_param)
  if(save.images) dev.off()
  # boxplot of thinning
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'GessThinBoxplot.pdf'))
  print(curve$gess_samplesThin)
  if(save.images) dev.off()
  
  # getting all the draws from the gess
  draws_list <- lapply(files[str_detect(files, paste0('.*',i,'.*', method,'_metrics'))], extract_draws)
  # plotting density of draws of gess
  if(save.images) pdf(file = paste0(file_path, 'curve',i, method,'Density.pdf'))
  facet_plot(draws_list)
  if(save.images) dev.off()
  # remove excess data
  rm(draws_list)
  gc()
  
  
  ### Transform
  method <- 'transform'
  
  # boxplot of x values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'TransformXBoxplot.pdf'))
  print(curve$transform_x)
  if(save.images) dev.off()
  # boxplot inv cdf values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'TransformInvCdfBoxplot.pdf'))
  print(curve$transform_inv_cdf)
  if(save.images) dev.off()
  # boxplot of parameters values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'TransformParamBoxplot.pdf'))
  print(curve$transform_param)
  if(save.images) dev.off()
  # boxplot of thining values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'TransformThinBoxplot.pdf'))
  print(curve$transform_samplesThin)
  if(save.images) dev.off()
  # scatter plot of kld
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'TransformKldPlot.pdf'))
  print(curve$transform_kld)
  if(save.images) dev.off()
  
  # gettig the draws from the transform slice sampler
  draws_list <- lapply(files[str_detect(files, paste0('.*',i,'.*', method,'_metrics'))], extract_draws)
  # plotin the density of the draws
  if(save.images) pdf(file = paste0(file_path, 'curve',i, method,'Density.pdf'))
  facet_plot(draws_list)
  if(save.images) dev.off()
  # removing excess data
  rm(draws_list)
  gc()
  
  
  ### Random Walk
  method <- 'rand_walk'
  
  # boxplot of x values
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'RandomWalkXBoxplot.pdf'))
  print(curve$rand_walk_x)
  if(save.images) dev.off()
  # boxplot of c values, sd of the envelope function
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'RandomWalkCBoxplot.pdf'))
  print(curve$rand_walk_c)
  if(save.images) dev.off()
  # boxplot of the combination of parameters
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'RandomWalkParamBoxplot.pdf'))
  print(curve$rand_walk_param)
  if(save.images) dev.off()
  # boxplot of thinning
  if(save.images) pdf(file = paste0(file_path, 'curve',i,'RandomWalkThinBoxplot.pdf'))
  print(curve$rand_walk_samplesThin)
  if(save.images) dev.off()
  
  # getting draws from the randwalk 
  draws_list <- lapply(files[str_detect(files, paste0('.*',i,'.*', method,'_metrics'))], extract_draws)
  # plotting the density of the draws
  if(save.images) pdf(file = paste0(file_path, 'curve',i, method,'Density.pdf'))
  facet_plot(draws_list)
  if(save.images) dev.off()
  # removing excess data
  rm(draws_list, curve)
  gc()

}

