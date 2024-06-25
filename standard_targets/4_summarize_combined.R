rm(list=ls())
library("tidyverse")

targets <- "all"
targets <- c("normal", "gamma", "igamma")
# targets <- c("normal", "gamma", "gammalog", "igamma", "igammalog")
targets <- c("gamma", "gammalog", "igamma", "igammalog")

rnd <- 2
dte <- 240229 # 10 parallel jobs
# dte <- 240301 # 20 parallel jobs
dte <- 240308 # 10 parallel jobs
dte <- 240309 # 10 parallel jobs, all targets randomized together
dte <- 240520 # 10 parallel jobs, all targets randomized together; qslice package

if (targets == "all") {
  dat <- read.csv(paste0("output/combined_round", rnd, "_all_", dte, ".csv"))
} else {
  datl <- list()
  for (tg in targets) {
    datl[[tg]] <- read.csv(paste0("output/combined_round", rnd, "_target", tg, "_", dte, ".csv"))
  }
  dat <- do.call(rbind, datl)
}

str(dat)
hist(dat$ks_pval)

## summarize Kolmogorov-Smirnov test statistics
dat$type_sub <- paste(dat$type, dat$subtype)
ks_rej <- dat %>% group_by(type_sub, target) %>% summarize(ks_rej = mean(ks_pval < 0.05), n = n())
ggplot(ks_rej, aes(x = jitter(ks_rej), y = type_sub, col = target)) + geom_point()

summary(ks_rej)


## summarize algorithm settings
algo_descrips <- dat %>% filter(!grepl("samples|MM", subtype)) %>%
	group_by(type_sub, target) %>% summarize(n = n(), t = length(unique(algo_descrip)))
print(algo_descrips, n = 50)
dat %>% filter(!grepl("samples|MM", subtype), target %in% c("normal", "gamma", "igamma")) %>%
	group_by(type_sub, target) %>% summarize(n = n(), t = unique(algo_descrip)) %>% print(., n = 50)

dat <- dat %>% mutate(algo = paste(type, subtype), algo_all = paste(type, algo_descrip))
dat$algo <- gsub(" NA", "", x = dat$algo)
unique(dat$algo)
unique(dat$algo_al)

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
  factor(., levels = c("Rand walk", "Slice", "Quantile slice", "IMH"))

dat$targetlab <- case_match(dat$target, "normal" ~ "Normal",
                            "gamma" ~ "Gamma",
                            "gammalog" ~ "Gamma-log",
                            "igamma" ~ "Inverse Gamma",
                            "igammalog" ~ "Inverse Gamma-log")

dat$target_base <- case_match(dat$target, "normal" ~ "Normal",
                              c("gamma", "gammalog") ~ "Gamma",
                              c("igamma", "igammalog") ~ "Inverse Gamma") %>%
  factor(., levels = c("Normal", "Gamma", "Inverse Gamma"))

dat$target_tx <- ifelse(grepl("log", dat$target), "Log transform", "Original") %>%
  as.factor()
dat$target_tx_alpha <- ifelse(dat$target_tx == "Log transform", 0.5, 1.0)


plt <- ggplot(dat %>% filter(target %in% c("normal", "gamma", "igamma")),
              aes(x = sampPsec / 1e3, y = algoF, fill = typeF), color = "gray") +
  # geom_boxplot() +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Effective samples per second\n(thousands)") + ylab("") + labs(color = "", fill = "") +
  facet_wrap(~ target_base, scales = "free_x")

plt

ggsave(plot = plt, width = 9, height = 5.5,
       filename = paste0("plots/ESPS_primary_round", rnd, "_", dte, ".pdf"))


plt <- ggplot(dat %>% filter(target %in% c("gamma", "igamma", "gammalog", "igammalog")),
              aes(x = sampPsec / 1e3, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  # geom_boxplot() +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Effective samples per second\n(thousands)") + ylab("") + labs(color = "", fill = "") +
  facet_wrap(~ target_base, scales = "fixed") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 9, height = 8,
       filename = paste0("plots/ESPS_gammas_w_tnx_round", rnd, "_", dte, ".pdf"))




