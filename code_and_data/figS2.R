###################################################

rm(list = ls(all = TRUE))

###################################################
## Basic
###################################################

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(CURRENT_WORKING_DIR)

load("kakunin/MCMC_mixture.Rdata")

mcmc_pars <- as.data.frame(MCMC_mixture$pars)

names(mcmc_pars)

# burn-in を除く場合：最初の20%を除外
burn <- floor(0.2 * nrow(mcmc_pars))
mcmc_use <- mcmc_pars[(burn + 1):nrow(mcmc_pars), ]

# original scale での相関
cor_raw <- cor(mcmc_use$lambda, mcmc_use$rr, method = "pearson", use = "complete.obs")

# log10 scale での相関
cor_log <- cor(log10(mcmc_use$lambda), log10(mcmc_use$rr), method = "pearson", use = "complete.obs")

cor_raw
cor_log

library(ggplot2)

df_lr <- data.frame(
  lambda = mcmc_use$lambda,
  r = mcmc_use$rr,
  log10_lambda = log10(mcmc_use$lambda),
  log10_r = log10(mcmc_use$rr)
)

cor_log <- cor(df_lr$log10_lambda, df_lr$log10_r, method = "pearson")

p <- ggplot(df_lr, aes(x = log10_lambda, y = log10_r)) +
  geom_point(alpha = 0.25, size = 0.7) +
  labs(
    x = expression(log[10](lambda)),
    y = expression(log[10](r)),
    title = paste0("Pearson's r = ", round(cor_log, 3))
  ) +
  theme_classic(base_size = 14)

ggsave("supplementary/Supplementary_Fig_S2_lambda_r_correlation.pdf", p, width = 5, height = 4)
ggsave("supplementary/Supplementary_Fig_S2_lambda_r_correlation.png", p, width = 5, height = 4, dpi = 300)
