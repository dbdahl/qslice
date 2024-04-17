rm(list=ls())
library("tidyverse")

target <- "normal"
target <- "gamma"
target <- "gammalog"
target <- "igamma"
target <- "igammalog"

rnd <- 1
dte <- 240229 # 10 parallel jobs
# dte <- 240301 # 20 parallel jobs
dte <- 240308 # 10 parallel jobs

dat <- read.csv(paste0("output/combined_round", rnd, "_target", target, "_", dte, ".csv"))

str(dat)
hist(dat$ks_pval)

dat <- dat %>% mutate(algo = paste(type, subtype), algo_all = paste(type, algo_descrip))
dat$algo <- gsub(" NA", "", x = dat$algo)
unique(dat$algo)
unique(dat$algo_al)

(medESPS_AUC <- median(dat$sampPsec[dat$algo == "Qslice AUC"]))

dat$algoF <- factor(dat$algo, levels = rev(c("rw", "imh AUC", "imh AUC_wide",
                                         "stepping", "gess", "latent",
                                         "Qslice MSW", "Qslice AUC",
                                         "Qslice MSW_samples", "Qslice AUC_samples",
                                         "Qslice AUC_wide",
                                         "Qslice Laplace_Cauchy", "Qslice MM_Cauchy")),
                    labels = rev(c("Random walk", "AUC", "AUC-diffuse",
                               "Step & shrink", "Gen. elliptical", "Latent",
                               "MSW", "AUC",
                               "MSW-samples", "AUC-samples",
                               "AUC-diffuse",
                               "Laplace-Cauchy", "MM-Cauchy"))
                    )

dat$typeF <- case_match(dat$type, c("gess", "latent", "stepping") ~ "Slice",
                        "imh" ~ "IMH", "Qslice" ~ "Quantile slice",
                        "rw" ~ "Rand walk") %>%
  factor(., levels = c("Rand walk", "IMH", "Slice", "Quantile slice"))

dat$targetlab <- case_match(dat$target, "normal" ~ "normal",
                            "gamma" ~ "gamma",
                            "gammalog" ~ "gamma-log",
                            "igamma" ~ "inverse gamma",
                            "igammalog" ~ "inverse gamma-log")

plt <- ggplot(dat, aes(x = sampPsec / 1e3, y = algoF, fill = typeF)) +
  geom_violin(draw_quantiles = 0.5) +
  # geom_boxplot() +
  ggtitle(target) + theme_bw() +
  xlab("Effective samples per second\n(thousands)") + ylab("") + labs(color = "", fill = "") # +
#  scale_x_continuous(sec.axis = sec_axis(~ . / medESPS_AUC, name = "ESPS relative to AUC median"))

plt

# ggsave(plot = plt, width = 6, height = 6,
#        filename = paste0("plots/ESPS_", target, "_round", rnd, "_", dte, ".pdf"))

ggplot(dat %>% filter(type == "rw"), aes(x = sampPsec, y = algo_all, fill = type)) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point() +
  ggtitle(target)

