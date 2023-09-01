# this is a summary of the results from the samplers with a gamma target
# author: Sam Johnson

# sourcing the files the different samplers
filesToSource <- Sys.glob('../../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)

library(tidyverse)

theme_set(
  theme_minimal() + 
    theme(legend.position = '',
          panel.grid = element_blank())
)


# files <- Sys.glob('data/*.rds')
# output <- sapply(files, FUN = readRDS)
targetTitle <- 'Gamma Target'

transform_metrics <- read_in_data('output/transform/*.csv')
stepping_out_metrics <- read_in_data('output/steppingout/*.csv')
latent_metrics <- read_in_data('output/latent/*.csv')
gess_metrics <- read_in_data('output/gess/*.csv')
rand_walk_metrics <- read_in_data('output/randwalk/*.csv')

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

pdf(file = 'images/transform_sampPsec.pdf')
transform_metrics %>% 
  ggplot(aes(y = unlist(t), x = sampPsec)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'psuedo target'
  )
dev.off()

pdf(file = 'images/transform_sampPsecAndnEval.pdf')
transform_metrics %>% 
  ggplot(aes(x = nEval, y = sampPsec, color = unlist(t))) + 
  geom_point(alpha = 0.5) + 
  geom_point(data = transform_metrics %>% 
               group_by(unlist(t)) %>% 
               summarise(nEval = mean(nEval), sampPsec = mean(sampPsec)) %>% 
               rename('t' = 'unlist(t)'),
             aes(x = nEval, y = sampPsec, color = t), size = 5, shape = 9)
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

pdf(file = 'images/steppingOut_sampPsec.pdf')
stepping_out_metrics %>% 
  ggplot(aes(y = factor(w), x = sampPsec)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'w'
  )
dev.off()

pdf(file = 'images/steppingOut_sampPsecAndnEval.pdf')
stepping_out_metrics %>% 
  dplyr::filter(w > 3) %>% 
  ggplot(aes(x = nEval, y = sampPsec, color = factor(w))) +
  geom_point(alpha = 0.5) +
  geom_point(data = stepping_out_metrics %>% 
               dplyr::filter(w > 3) %>% 
               group_by(w) %>% 
               summarise(sampPsec = mean(sampPsec), nEval = mean(nEval)),
             aes(x = nEval, y = sampPsec, color = factor(w)), size = 5, shape = 9)
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

pdf(file = 'images/latent_sampPsec.pdf')
latent_metrics %>% 
  ggplot(aes(y = factor(rate), x = sampPsec)) +
  geom_boxplot() + 
  labs(
    title = targetTitle,
    y = 'rate'
  )
dev.off()

pdf(file = 'images/latent_sampPsecAndnEval.pdf')
latent_metrics %>% 
  ggplot(aes(x = nEval, y = sampPsec, color = factor(rate))) +
  geom_point(alpha = 0.5) +
  geom_point(data = latent_metrics %>% 
               group_by(rate) %>% 
               summarise(sampPsec = mean(sampPsec), nEval = mean(nEval)),
             aes(x = nEval, y = sampPsec, color = factor(rate)), size = 5, shape = 9)
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

pdf(file = 'images/gess_sampPsec.pdf')
gess_metrics %>% 
  ggplot(aes(y = paste0(mu,',',sigma,',',df), x = sampPsec)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'mu,sigma,df'
  )
dev.off()

pdf(file = 'images/gess_sampPsecAndnEval.pdf')
gess_metrics %>%
  mutate(par = paste0(mu,',',sigma,',',df)) %>% 
  ggplot(aes(x = nEval, y = sampPsec, color = factor(par))) +
  geom_point(alpha = 0.5) +
  geom_point(data = gess_metrics %>% 
               mutate(par = paste0(mu,',',sigma,',',df)) %>% 
               group_by(par) %>% 
               summarise(sampPsec = mean(sampPsec), nEval = mean(nEval)),
             aes(x = nEval, y = sampPsec, color = factor(par)), size = 5, shape = 9)
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

pdf(file = 'images/randWalk_sampPsec.pdf')
rand_walk_metrics %>% 
  ggplot(aes(y = factor(c), x = sampPsec)) +
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = 'c'
  )
dev.off()

pdf(file = 'images/randWalk_AccptanceRate_sampPsec.pdf')
rand_walk_metrics |> 
  ggplot(aes(x = acceptanceRate, y = sampPsec, color = factor(c))) + 
  geom_point() + 
  theme(
    legend.position = 'bottom'
  ) +
  geom_vline(xintercept = 0.43, linetype = 'dashed', color = 'red') + 
  annotate("text", x = .41, y = 5000, label = "43%", color = 'red')
dev.off()

# combining them all
totalsampPsec <- combine_func(stepping_out_metrics = stepping_out_metrics, 
                               latent_metrics = latent_metrics,
                               gess_metrics = gess_metrics,
                               rand_walk_metrics = rand_walk_metrics,
                               transform_metrics = transform_metrics,
                               metric = 'sampPsec') |> 
  rowwise() |> 
  mutate(waterPenatly = case_when(
    grepl('*c2:0.01 [A-Z]*',par) ~ 0.01,
    grepl('*c2:0.1 [A-Z]*',par) ~ 0.1,
    grepl('*c2:0.5 [A-Z]*',par) ~ 0.5,
    grepl('*c2:0.66 [A-Z]*',par) ~ 0.66,
    grepl('*c2:0.75 [A-Z]*',par) ~ 0.75,
    grepl('*c2:0 [A-Z]*',par) ~ 0,
    grepl('*c2:1.33 [A-Z]*',par) ~ 1.33,
    grepl('*c2:1.5 [A-Z]*',par) ~ 1.5,
    grepl('*c2:2 [A-Z]*',par) ~ 2,
    grepl('*c2:100 [A-Z]*',par) ~ 100,
    grepl('*c2:10 [A-Z]*',par) ~ 10,
    grepl('*c2:1 [A-Z]*',par) ~ 1,
    TRUE ~ 99
  ),
  optimMethod = case_when(
    grepl('*OSAUC', par) ~ 'OSAUC',
    grepl('*OAUC', par) ~ 'OAUC',
    grepl('*OS', par) ~ 'OS',
    grepl('*O', par) ~ 'O',
    grepl('*Man|*Auto|*Laplace|*MM', par) ~ 'Man',
    grepl('*w=*|*c=*|*rate=*|*mu=*', par) ~ 'Comp',
    TRUE ~ 'Other'
  ))

Comp <- totalsampPsec$par[!grepl('*loc = *', totalsampPsec$par)] |> unique()
Transform <- totalsampPsec$par[grepl('.*Man|.*MM|.*Auto|.*Laplace', totalsampPsec$par)] |> unique()
OSAUC <- totalsampPsec$par[grepl('*OSAUC*', totalsampPsec$par)] |> unique()
OAUC <- totalsampPsec$par[grepl('*OAUC', totalsampPsec$par)] |> unique()
O <- totalsampPsec$par[grepl('.*O$', totalsampPsec$par)] |> unique()
OS <- totalsampPsec$par[grepl('.*OS$', totalsampPsec$par)] |> unique()

totalsampPsec$par <- factor(totalsampPsec$par, levels = c(Comp, Transform, OSAUC, OAUC, O, OS))

# mutate(par = forcats::fct_reorder2(par,optimMethod,waterPenatly))
# mutate(par = fct_reorder(par,optimMethod))

pdf(file = 'images/totalsampPsec.pdf')
totalsampPsec %>% 
  dplyr::filter(optimMethod != 'Transform') |> 
  ggplot(aes(y = par, x = sampPsec, color = method, fill = optimMethod)) + 
  geom_boxplot() +
  labs(
    title = targetTitle
  )
dev.off()

## simplified
pdf(file = 'images/gammatotalsampPsecSimplfied.pdf')
totalsampPsec |> 
  dplyr::filter(optimMethod %in% c('OSAUC','OAUC','Man','Comp','Transform')) |> 
  dplyr::filter(!grepl('*Auto|*Man',par)) |> 
  dplyr::filter(waterPenatly == 0 || waterPenatly == 99) |> 
  mutate(par = case_when(
    grepl('*mu=*', par) ~ 'GESS',
    grepl('*rate=*', par) ~ 'Latent',
    grepl('*MM', par) ~ 'Moment Matching',
    grepl('*w=*', par) ~ 'Stepping Out',
    grepl('*Laplace*', par) ~ 'Laplace',
    grepl('*OSAUC*', par) ~ 'AUC Samples',
    grepl('*OAUC*', par) ~ 'AUC',
    grepl('*c=*', par) ~ 'Random Walk'
  )) |>
  mutate(par = ordered(par, levels = c('Stepping Out','GESS','Latent','Random Walk','Laplace','Moment Matching','AUC Samples','AUC'))) |> 
  ggplot(aes(y = par, x = sampPsec, fill = method)) + 
  geom_boxplot() +
  labs(
    title = targetTitle,
    y = '',
    x = 'Effective Samples Per CPU Second'
  ) +
  scale_x_continuous(labels = scales::number_format(scale = 1e-3, suffix = 'K'))
dev.off()

