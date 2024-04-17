rm(list=ls()); dev.off()
library("tidyverse")

targets <- "all"

dte <- 240330
dte <- 240405 # extra computation in full conditional


if (targets == "all") {
  dat <- read.csv(paste0("output/combined_all_", dte, ".csv"))
} else {
  datl <- list()
  for (tg in targets) {
    datl[[tg]] <- read.csv(paste0("output/combined_target", tg, "_", dte, ".csv"))
  }
  dat <- do.call(rbind, datl)
}

str(dat)


## summarize algorithm settings
dat <- dat %>% mutate(algo = paste(type, subtype))
dat$algo <- gsub(" NA", "", x = dat$algo)
unique(dat$algo)

dat$algoF <- factor(dat$algo, levels = rev(c("rw", "stepping", "latent",
                                             "imh AUC_samples", "imh Laplace_analytic", "imh Laplace_analytic_wide",
                                             "gess AUC_samples", "gess Laplace_analytic", "gess Laplace_analytic_wide",
                                             "Qslice AUC_samples", "Qslice Laplace_analytic", "Qslice Laplace_analytic_wide")),
                    labels = rev(c("Random walk", "Step & shrink", "Latent",
                                   "Pseudo: AUC", "Pseudo: Laplace", "Pseudo: Laplace (wide)",
                                   "Pseudo: AUC", "Pseudo: Laplace", "Pseudo: Laplace (wide)",
                                   "Pseudo: AUC", "Pseudo: Laplace", "Pseudo: Laplace (wide)"))
)

dat$typeF <- case_match(dat$type, c("gess", "latent", "stepping") ~ "Slice",
                        "imh" ~ "IMH", "Qslice" ~ "Quantile slice",
                        "rw" ~ "Rand walk") %>%
  factor(., levels = c("Rand walk", "Slice", "Quantile slice", "IMH"))

dat$targetlab <- case_match(dat$target, "hyper-g" ~ "h-g",
                            "hyper-g-log" ~ "h-g log")

dat$target_tx <- ifelse(grepl("log", dat$target), "Log transform", "Original") %>%
  as.factor()
dat$target_tx_alpha <- ifelse(dat$target_tx == "Log transform", 0.5, 1.0)


plt <- ggplot(dat %>% filter(target %in% c("hyper-g")),
              aes(x = sampPsec / 1e3, y = algoF, fill = typeF), color = "gray") +
  geom_violin(draw_quantiles = 0.5 , scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Effective samples per second\n(thousands)") + ylab("") + labs(color = "", fill = "")

plt

ggsave(plot = plt, width = 6, height = 6,
       filename = paste0("plots/ESPS_primary_", dte, ".pdf"))


plt <- ggplot(dat %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = sampPsec / 1e3, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  # geom_boxplot() +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Effective samples per second\n(thousands)") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 6, height = 6,
       filename = paste0("plots/ESPS_w_tnx_", dte, ".pdf"))



## diagnostic

plt <- ggplot(dat %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = g_mn, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Posterior mean: g") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/pm_g", dte, ".pdf"))

plt <- ggplot(dat %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = g_sd, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Posterior SD: g") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

plt <- ggplot(dat %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = g_SE, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") + xlim(c(0, 0.2)) +
  xlab("MCMC std. error: g") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/SE_g", dte, ".pdf"))

plt <- ggplot(dat %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = nEval / n_iter, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Target evaluations per iteration") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

dat %>% group_by(type, subtype) %>% summarize(mn_eval = mean(nEval/n_iter), sd_eval = sd(nEval/n_iter))

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/nEval_", dte, ".pdf"))

plt <- ggplot(dat %>% filter(target %in% c("hyper-g", "hyper-g-log"), type == "Qslice"),
              aes(x = auc, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Transformed pseudo-target: AUC") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 6, height = 5,
       filename = paste0("plots/AUC_", dte, ".pdf"))
