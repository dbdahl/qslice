rm(list=ls()); dev.off()
library("tidyverse")

# targets <- "all"
targets <- "db"
targets <- "db40"
# targets <- c("db40", "db")

# dte <- 260518
# dte <- 260520 # includes samples_reg
dte <- 260521 # includes samples_reg with pseudo_opt on residuals

reference_algo <- "rw"
reference_tnx <- "tau2_marg-log" # $ at end means no -log appended

if (length(targets) > 0) {
  datl <- list()
  for (tg in targets) {
    datl[[tg]] <- read.csv(paste0("output/combined_", tg, "_", dte, ".csv"))
    med_sampPsec_ref <- median(filter(datl[[tg]], type == reference_algo, grepl(reference_tnx, target))$sampPsec)
    med_sampPiter_ref <- median( with(filter(datl[[tg]], type == reference_algo, grepl(reference_tnx, target)), ESS / n_iter) )
    datl[[tg]]$sampPsec_rel <- datl[[tg]]$sampPsec / med_sampPsec_ref
    datl[[tg]]$sampPiter_rel <- (datl[[tg]]$ESS / datl[[tg]]$n_iter)  / med_sampPiter_ref
  }
  dat <- do.call(rbind, datl)
} else if (targets == "all") {
  dat <- read.csv(paste0("output/combined_all_", dte, ".csv"))
}

str(dat)


## hack to create a placeholder for run combinations that did not occur

for (i in 1:2) {
  dat[nrow(dat) + 1, ] <- NA
  dat[nrow(dat), c("target", "type", "subtype")] <- c("db40_tau2_marg", "imh", "samples_reg")
  dat[nrow(dat), c("sampPsec_rel")] <- dat[nrow(dat), c("sampPiter_rel")] <- 1.0
}
for (i in 1:2) {
  dat[nrow(dat) + 1, ] <- NA
  dat[nrow(dat), c("target", "type", "subtype")] <- c("db40_tau2_marg", "gess", "samples_reg")
  dat[nrow(dat), c("sampPsec_rel")] <- dat[nrow(dat), c("sampPiter_rel")] <- 1.0
}
for (i in 1:2) {
  dat[nrow(dat) + 1, ] <- NA
  dat[nrow(dat), c("target", "type", "subtype")] <- c("db40_tau2_marg", "Qslice", "samples_reg")
  dat[nrow(dat), c("sampPsec_rel")] <- dat[nrow(dat), c("sampPiter_rel")] <- 1.0
}
for (i in 1:2) {
  dat[nrow(dat) + 1, ] <- NA
  dat[nrow(dat), c("target", "type", "subtype")] <- c("db_tau2_marg", "imh", "samples_reg")
  dat[nrow(dat), c("sampPsec_rel")] <- dat[nrow(dat), c("sampPiter_rel")] <- 1.0
}
for (i in 1:2) {
  dat[nrow(dat) + 1, ] <- NA
  dat[nrow(dat), c("target", "type", "subtype")] <- c("db_tau2_marg", "gess", "samples_reg")
  dat[nrow(dat), c("sampPsec_rel")] <- dat[nrow(dat), c("sampPiter_rel")] <- 1.0
}
for (i in 1:2) {
  dat[nrow(dat) + 1, ] <- NA
  dat[nrow(dat), c("target", "type", "subtype")] <- c("db_tau2_marg", "Qslice", "samples_reg")
  dat[nrow(dat), c("sampPsec_rel")] <- dat[nrow(dat), c("sampPiter_rel")] <- 1.0
}


tail(dat, n = 10)



dat$n <- ifelse(grepl("40", dat$target), 40, 442)


## summarize algorithm settings
dat <- dat %>% mutate(algo = paste(type, subtype))
dat$algo <- gsub(" NA", "", x = dat$algo)
unique(dat$algo)

dat$algoF <- factor(dat$algo, levels = rev(c("rw", "stepping", "latent",
                                             "imh AUC_samples", "imh MSW_samples", "imh samples_reg",
                                             "gess AUC_samples", "gess MSW_samples", "gess samples_reg",
                                             "Qslice AUC_samples", "Qslice MSW_samples", "Qslice samples_reg")),
                    labels = rev(c("Random walk", "Steping out & shrinkage", "Latent slice",
                                   "Independence M-H: AUC", "Independence M-H: MSW", "Independence M-H: reg",
                                   "Generalized elliptical slice: AUC", "Generalized elliptical slice: MSW", "Generalized elliptical slice: reg",
                                   "Quantile slice: AUC", "Quantile slice: MSW", "Quantile slice: reg"))
)


dat$typeF <- case_match(dat$type, c("gess", "latent", "stepping") ~ "Slice",
                        "imh" ~ "IMH", "Qslice" ~ "Quantile slice",
                        "rw" ~ "Rand walk") %>%
  factor(., levels = c("Rand walk", "Slice", "Quantile slice", "IMH"))

dat$target_tx <- ifelse(grepl("log", dat$target), "Log transform", "Original") %>%
  as.factor()
dat$target_tx_alpha <- ifelse(dat$target_tx == "Log transform", 0.5, 1.0)


## evaluations per iteration
eval_tab <- dat %>% group_by(target, type, subtype) %>% summarize(eval_per_iter_mean = mean(nEval / n_iter),
                                                      eval_per_iter_sd = sd(nEval / n_iter))

print(eval_tab, n = 30)

tail(dat, n = 10)




dat_plt <- dat %>% filter(n %in% c(40))
dat_plt <- dat %>% filter(n %in% c(442))

## tuning parameters
ggplot(dat %>% filter(grepl("tau2_marg$", target), algo == "latent"),
       aes(x = tuneParam, y = algoF, fill = typeF)) +
  geom_violin() +  theme_bw() + theme(legend.position = "none") +
  xlab("Tuning parameter value") + ylab("") + labs(color = "", fill = "")

ggplot(dat %>% filter(grepl("tau2_marg-log", target), algo == "latent"),
       aes(x = tuneParam, y = algoF, fill = typeF)) +
  geom_violin() +  theme_bw() + theme(legend.position = "none") +
  xlab("Tuning parameter value") + ylab("") + labs(color = "", fill = "")

ggplot(dat %>% filter(grepl("tau2_marg-log", target), algo == "Qslice samples_reg"),
       aes(x = tuneParam, y = algoF, fill = typeF)) +
  geom_violin() +  theme_bw() + theme(legend.position = "none") +
  xlab("Tuning parameter value") + ylab("") + labs(color = "", fill = "")



## timing results
plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg$", target)),
              aes(x = sampPsec / median(sampPsec[which(algo == "stepping")]), y = algoF, fill = typeF), color = "gray") +
  geom_violin(draw_quantiles = 0.5 , scale = "width") +
  theme_bw() + theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples per second") + ylab("") + labs(color = "", fill = "")

plt

ggsave(plot = plt, width = 6, height = 6,
       filename = paste0("plots/ESPS_primary_", dte, ".pdf"))


plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = sampPsec / median(sampPsec[which(algo == "stepping" & grepl("tau2_marg$", target))]),
                  y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  # geom_boxplot() +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  expand_limits(x = 0) +
  xlab("Effective samples per second") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 6, height = 6,
       filename = paste0("plots/ESPS_w_tnx_", dte, ".pdf"))






## add evals per iter on paneled plot

n_select <- 40
# n_select <- 442
# n_select <- c(40, 442)

dat_plt_n <- dat %>% filter(n %in% n_select) %>%
  mutate(n_lab = factor(n,
                        levels = c(40, 442),
                        labels = c("n = 40, p = 64", "n = 442, p = 64")))


format_mean_sd <- function(mn, std) {
  # Round std for conditional logic
  std_fmt <- ifelse(std < 0.005 & std != 0, 0.01, std)

  # Choose format string
  ifelse(
    std == 0,
    sprintf("%.0f ± %.0f", mn, std),
    sprintf("%.1f ± %.2f", mn, std_fmt)
  )
}

stat_labels <- dat_plt_n %>%
  filter(grepl("tau2_marg", target)) %>%
  group_by(algoF, target_tx, n_lab) %>%
  summarize(mEval = mean(nEval / n_iter),
            sdEval = sd(nEval / n_iter),
            maxSampPsec = max(sampPsec_rel),
            .groups = "drop") %>%
  mutate(algoF_num = as.numeric(factor(algoF, levels = sort(unique(dat_plt_n$algoF)))),
         dodge_offset = ifelse(target_tx == "Original", 0.2, -0.2),
         y = algoF_num + dodge_offset,
         x = Inf,
         label = format_mean_sd(mEval, sdEval)
  )

# x_pad <- dat_plt_n %>% # hack to start the horizontal axis at 0
#   group_by(n_lab) %>%
#   slice(1) %>%  # grab one row per facet
#   mutate(sampPsec = 0)  # force x = 0
#
# dat_plt_padded <- bind_rows(dat_plt_n, x_pad)

plt <- ggplot(dat_plt_n %>% filter(grepl("tau2_marg", target)),
              aes(x = sampPsec_rel,
                  y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha,
                  group = interaction(algoF, target_tx, drop = FALSE))) +
  expand_limits(x = 0) +
  geom_vline(xintercept = 1.0, col = "black") +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples per second") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black")) +
  geom_text(data = stat_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3)

plt_n <- plt +
  facet_grid(cols = vars(n_lab), scales = "fixed") +
  expand_limits(x = 0) +
  coord_cartesian(clip = "off") + theme(panel.spacing = unit(3.5, "lines")) +
  theme(plot.margin = unit(c(5.5, 50, 5.5, 0), "pt")) +
  theme(panel.border = element_blank()) +
  # geom_blank(data = x_pad, aes(x = sampPsec_rel)) +
  theme(strip.background = element_rect(fill = "gray90", color = NA))
plt_n

ggsave(plot = plt_n, width = 8, height = 4.5,
       filename = paste0("plots/hs_ESPS_w_tnx_n_", dte, ".pdf"))



(xlim <- range(c(range(dat_plt_n$sampPsec_rel), range(dat_plt_n$sampPiter_rel))))

plt_esps <- ggplot(dat_plt_n %>% filter(grepl("tau2_marg", target)),
                   aes(x = sampPsec_rel,
                       y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha,
                       group = interaction(algoF, target_tx, drop = FALSE))) +
  geom_vline(xintercept = 1.0, col = "gray20") +
  geom_blank() +
  geom_violin(draw_quantiles = 0.5, scale = "width", position = "dodge", drop = FALSE) +
  theme_bw() + theme(legend.position = "none") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt")) +
  theme(panel.border = element_blank()) +
  scale_alpha(guide = "none") +
  xlab("Effective samples per second") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "gray20"), drop = FALSE)
# plt_esps

plt_espit0 <- ggplot(dat_plt_n %>% filter(grepl("tau2_marg", target)),
                   aes(x = sampPiter_rel,
                       y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha,
                       group = interaction(algoF, target_tx, drop = FALSE))) +
  geom_vline(xintercept = 1.0, col = "gray20") +
  geom_blank() +
  geom_violin(draw_quantiles = 0.5, scale = "width", position = "dodge", drop = FALSE) +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_alpha(guide = "none") +
  xlab("Effective samples per iteration (1/IAT)") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "gray20")) +
  geom_text(data = stat_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = "outward", size = 3.25)

plt_espit <- plt_espit0 +
  coord_cartesian(clip = "off") + theme(panel.spacing = unit(3.5, "lines")) +
  theme(plot.margin = unit(c(5.5, 55, 5.5, 0), "pt")) +
  theme(panel.border = element_blank()) +
  # geom_blank(data = x_pad, aes(x = sampPsec_rel)) +
  theme(strip.background = element_rect(fill = "gray90", color = NA))

# plt_espit

library("patchwork")

plt_panel <- plt_esps + plt_espit & scale_x_continuous(limits = xlim, breaks = c(0.5, 1.0, 1.5, 2.0, 2.5))

ggsave(plot = plt_panel, width = 8, height = 4.5,
       filename = paste0("plots/hs_esps_espit_tnx_n", n_select, ".pdf"))






dat_plt <- dat %>% filter(n %in% c(40))
dat_plt <- dat %>% filter(n %in% c(442))
dat_plt <- dat %>% filter(n %in% c(40, 442))



## diagnostic

plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = tau2_mn, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") + # xlim(c(0.0010, 0.004)) +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Posterior mean: tau2") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/pm_tau2", dte, ".pdf"))

plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = tau2_sd, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") + # xlim(c(0, 0.01)) +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Posterior SD: tau2") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = tau2_SE, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") + # xlim(c(0, 0.0003)) +
  xlab("MCMC std. error: tau2") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt


plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = ltau2_SE, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") + # xlim(c(0, 0.0003)) +
  xlab("MCMC std. error: log(tau2)") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/SE_logtau2", dte, ".pdf"))

plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = sampPiter_rel, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") + # xlim(c(0, 0.0003)) +
  xlab("Effective samples per iteration (1 / IAT), relative") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/espit_logtau2", dte, ".pdf"))


plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target)),
              aes(x = nEval / n_iter, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Target evaluations per iteration") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

dat_plt %>% group_by(type, subtype) %>% summarize(mn_eval = mean(nEval/n_iter), sd_eval = sd(nEval/n_iter))

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/nEval_", dte, ".pdf"))

plt <- ggplot(dat_plt %>% filter(grepl("tau2_marg", target), type == "Qslice"),
              aes(x = auc, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Transformed pseudo-target: AUC") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 6, height = 5,
       filename = paste0("plots/AUC_", dte, ".pdf"))
