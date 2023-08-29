## A script to summarize the gprior results
# Author: Sam Johnson

filesToSource <- Sys.glob('../utilityFunctions/*.R')[-13]
sapply(filesToSource, source)
library(magrittr)
library(tidyverse)
library(forcats)
source('formatingFunctions.R')

theme_set(
  theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = '',
          plot.title = element_text(hjust = 0.5, size = 30),
          axis.title.x = element_text(size = 25), 
          axis.text = element_text(size = 20),
          panel.grid.major = element_blank()
          # axis.line.y.left = element_line(color = 'black'),
          # axis.line.x.bottom = element_line(color = 'black'),
          # axis.ticks = element_line(color = 'black'),
          # axis.ticks.length=unit(.25, "cm"),
          # plot.title = element_text(hjust = 0.5)
          
    )
)

# reading in all the files
# files <- (Sys.glob("data/*.rds"))
# data <- sapply(files, readRDS, simplify = FALSE)
# results <- lapply(data, resultsFunc)

files <- Sys.glob("output/*/*.csv")
data <- sapply(files, read.csv, simplify = FALSE)
results <- do.call(rbind, data)  |> 
  mutate(
    sampler = case_when(
      grepl('*Samples*|*auto*|*laplace*',t) ~ 'Transform',
      grepl('*latent*',t) ~ 'Latent',
      grepl('*gess*',t) ~ 'GESS',
      grepl('*steppingout*',t) ~ 'Stepping Out',
      grepl('*independence*',t) ~ 'Independence',
      grepl('*randwalk*',t) ~ 'Random Walk',
      TRUE ~ 'Other'
    ),
    pseudoType = case_when(
      grepl('*optimSamplesAUC*',t) ~ 'OSAUC',
      grepl('*optimSample*',t) ~ 'OS',
      grepl('*auto*',t) ~ 'Auto',
      grepl('*laplace*',t) ~ 'Laplace',
      TRUE ~ 'Other'
    ))

pdf(file = 'images/boxplot.pdf')
results |> 
  ggplot(aes(x = sampPsec, y = sampler, fill = t)) +
  geom_boxplot() +
  theme(
    legend.position = 'bottom'
  ) +
  labs(
    title = 'Zellner G Prior',
    y = '',
    x = 'Effective Samples Per CPU Second',
    fill = 'Pseudo Target Method'
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.title.x = element_text(size = 20)
  )
dev.off()

pdf(file = 'images/bestboxplot.pdf')
results |> 
  dplyr::filter(t %in% c('optimSamplesAUC','latent','steppingout','gess','independence','randwalk')) |> 
  ggplot(aes(x = sampPsec, y = sampler, fill = sampler)) +
  geom_boxplot() +
  labs(
    title = 'Zellner G Prior',
    y = '',
    x = 'Effective Samples Per CPU Second'
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.title.x = element_text(size = 20)
  )
dev.off()
