###################################################
## Basic
###################################################

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)

rm(list = ls(all = TRUE))

###################################################
## Packages
###################################################

library(dplyr)
library(tidyr)
library(FME)
library(ggplot2)
library(ggpmisc)
library(DescTools)
library(HDInterval)
library(patchwork)
library(deSolve)
library(coda)

###################################################
## Basic configuration
###################################################

mcmc_iter <- 200000

Tmin <- 0.0
Tmax <- 60.0
step_size <- 0.1
stime <- seq(Tmin, Tmax, step_size)

lbt <- 0.001
ubt <- 1000.0

pars <- c(
  x0 = 1.829e5,
  y0 = 1.0,
  dx = 0.193,
  epsilon = 0.999,
  lambda = 1.115e14,
  f = 0.0877,
  alpha = 1,
  rho = 0.1,
  rr = 1.4e-7,
  tme30 = 30
)

fit_datatypes <- c(
  "condition1_ccc", "condition1_dna",
  "condition2_ccc", "condition2_dna",
  "condition3_ccc", "condition3_dna",
  "condition4_ccc", "condition4_dna"
)

validation_datatypes <- c(
  "ccc5", "dna5",
  "ccc6", "dna6",
  "ccc7", "dna7",
  "ccc8", "dna8"
)

Tmin_plot <- 0
Tmax_plot <- Tmax

###################################################
## Plot settings
###################################################

pltsetb <- theme(legend.position = "none") +
  theme(plot.background = element_blank()) +
  theme(text = element_text("Helvetica")) +
  theme(axis.text = element_text(colour = "black", size = 15)) +
  theme(axis.ticks = element_line(colour = "black", linewidth = 1)) +
  theme(axis.ticks.length = unit(0.2, "cm")) +
  theme(axis.line.x = element_line(colour = "black", linewidth = 1)) +
  theme(axis.line.y = element_line(colour = "black", linewidth = 1)) +
  theme(axis.title = element_text(colour = "black", size = 15)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

caka      <- "#CC3366"
cakao     <- "#CC336650"
cao       <- "#3366CC"
caoo      <- "#3366CC50"
cmurasaki <- "#6600FF"

###################################################
## Utilities
###################################################

read_obs_file <- function(path, stime) {
  dsf <- read.csv(path, sep = "\t", header = FALSE)
  colnames(dsf) <- c("time", "value")
  ids <- vapply(dsf$time, function(x) which.min(abs(stime - x)), integer(1))
  cbind(ids = ids, value = dsf$value)
}

read_log10_file <- function(path, sep = "\t") {
  dsf <- read.csv(path, sep = sep, header = FALSE)
  data.frame(
    time = dsf[, 1],
    value = log10(pmax(dsf[, 2], 1e-12))
  )
}

compile_final_dll <- function() {
  dll_name <- paste0("final", .Platform$dynlib.ext)
  if ("final" %in% names(getLoadedDLLs())) {
    dyn.unload(dll_name)
  }
  suppressWarnings(file.remove("final.o", dll_name))
  system("R CMD SHLIB final.c")
  dyn.load(dll_name)
}

set_constr <- function(low_vec, up_vec) {
  n_vec <- length(low_vec)
  stopifnot(length(up_vec) == n_vec)

  AL <- diag(n_vec)
  AU <- -diag(n_vec)
  Amat <- rbind(AL, AU)
  bvec <- c(low_vec, -up_vec)

  list(
    low.val = low_vec,
    up.val = up_vec,
    Amat = Amat,
    bvec = bvec
  )
}

make_sim_df <- function(fitted, state_name, type_name) {
  data.frame(
    time = fitted$time,
    value = log10(pmax(fitted[[state_name]], 1e-12)),
    type = type_name
  )
}

make_ci_df <- function(mat, times) {
  mean_vec <- log10(pmax(apply(mat, 2, mean), 1e-12))
  low_vec  <- log10(pmax(apply(mat, 2, function(x) quantile(x, 0.025)), 1e-12))
  high_vec <- log10(pmax(apply(mat, 2, function(x) quantile(x, 0.975)), 1e-12))

  list(
    ribbon = data.frame(x = times, mean = mean_vec, ymin = low_vec, ymax = high_vec),
    low    = data.frame(x = times, y = low_vec),
    mean   = data.frame(x = times, y = mean_vec),
    high   = data.frame(x = times, y = high_vec)
  )
}

make_check_df <- function(line_df, obs_df) {
  line_df %>%
    mutate(x = round(x, 1)) %>%
    rename(time = x, estimated = y) %>%
    left_join(obs_df, by = "time") %>%
    rename(observed = value) %>%
    na.omit()
}

mcmc_plot_func <- function(c.list, d.list, ccc, dna) {
  ggplot() +
    geom_ribbon(data = c.list$ribbon, aes(x = x, ymin = ymin, ymax = ymax), fill = caoo) +
    geom_line(data = c.list$low,  aes(x = x, y = y), linewidth = 1, color = caoo) +
    geom_line(data = c.list$mean, aes(x = x, y = y), linewidth = 1, color = cao) +
    geom_line(data = c.list$high, aes(x = x, y = y), linewidth = 1, color = caoo) +
    geom_ribbon(data = d.list$ribbon, aes(x = x, ymin = ymin, ymax = ymax), fill = cakao) +
    geom_line(data = d.list$low,  aes(x = x, y = y), linewidth = 1, color = cakao) +
    geom_line(data = d.list$mean, aes(x = x, y = y), linewidth = 1, color = caka) +
    geom_line(data = d.list$high, aes(x = x, y = y), linewidth = 1, color = cakao) +
    geom_point(data = ccc, aes(x = time, y = value), colour = cao, size = 4) +
    geom_point(data = dna, aes(x = time, y = value), colour = caka, size = 4) +
    labs(x = "Days post infection", y = "ccc / HBV DNA (copies/well)") +
    scale_x_continuous(
      limits = c(Tmin_plot, Tmax_plot),
      breaks = seq(0, 60, by = 10)
    ) +
    scale_y_continuous(
      limits = c(0, 12),
      breaks = seq(0, 12, 2),
      labels = c(expression(10^0), expression(10^2), expression(10^4),
                 expression(10^6), expression(10^8), expression(10^10),
                 expression(10^12))
    ) +
    pltsetb
}

fig_plot_func <- function(sccc, sdna, dccc, ddna) {
  ggplot() +
    geom_path(data = sccc, aes(x = time, y = value), color = cao, linewidth = 1) +
    geom_path(data = sdna, aes(x = time, y = value), color = caka, linewidth = 1) +
    geom_point(data = dccc, aes(x = time, y = value), color = cao, size = 6) +
    geom_point(data = ddna, aes(x = time, y = value), color = caka, size = 6) +
    labs(x = "days", y = "ccc/HBV DNA (copies/well)") +
    scale_x_continuous(breaks = seq(0, 60, by = 10)) +
    scale_y_continuous(limits = c(0, 12)) +
    pltsetb
}

make_corr_plot <- function(df, point_color, outfile) {
  tmp_ccc <- CCC(df$estimated, df$observed, ci = "z-transform", conf.level = 0.95)
  ccc_val <- tmp_ccc$rho.c[, "est"]
  ccc_ci  <- tmp_ccc$rho.c[, c("lwr.ci", "upr.ci")]

  plt <- ggplot(df, aes(x = observed, y = estimated)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
    geom_point(size = 2, alpha = 0.9, color = point_color) +
    stat_correlation(use_label(c("r", "P")), method = "pearson", size = 4.5) +
    annotate(
      "text",
      x = 4.63,
      y = 9,
      label = sprintf(
        "CCC = %.3f (95%% CI: %.3f–%.3f)",
        ccc_val, ccc_ci[1], ccc_ci[2]
      ),
      hjust = 0,
      size = 4.5
    ) +
    labs(y = "Estimated value", x = "Observed value") +
    scale_y_continuous(
      limits = c(4.6, 9.4),
      breaks = seq(5, 9, 1),
      labels = c(expression(10^5), expression(10^6), expression(10^7),
                 expression(10^8), expression(10^9))
    ) +
    scale_x_continuous(
      limits = c(4.6, 9.4),
      breaks = seq(5, 9, 1),
      labels = c(expression(10^5), expression(10^6), expression(10^7),
                 expression(10^8), expression(10^9))
    ) +
    theme(axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)) +
    pltsetb

  ggsave(outfile, plt, width = 5, height = 5)
  plt
}

###################################################
## Data loading
###################################################

obs_data <- lapply(
  fit_datatypes,
  function(x) read_obs_file(file.path("data", paste0(x, ".txt")), stime)
)

fit_data_list <- lapply(
  fit_datatypes,
  function(x) read_log10_file(file.path("data", paste0(x, ".txt")), sep = "\t")
)
names(fit_data_list) <- fit_datatypes

validation_paths <- file.path("data_kakunin", paste0(validation_datatypes, ".csv"))
missing_files <- validation_paths[!file.exists(validation_paths)]
if (length(missing_files) > 0) {
  stop("Missing validation files: ", paste(missing_files, collapse = ", "))
}

validation_data_list <- lapply(
  validation_paths,
  function(path) read_log10_file(path, sep = ",")
)
names(validation_data_list) <- validation_datatypes

###################################################
## Compile and load DLL
###################################################

compile_final_dll()

###################################################
## ODE models
###################################################

ODEs <- function(pars) {
  x00 <- as.numeric(pars["x0"])
  y00 <- as.numeric(pars["y0"])

  rhs <- c(
    x   = x00, y   = y00,
    x1s = x00, y1s = y00,
    x2s = x00, y2s = y00,
    x3s = x00, y3s = y00
  )

  times <- seq(Tmin, Tmax, step_size)

  out <- ode(
    y = rhs,
    parms = pars,
    times = times,
    func = "derivs",
    initfunc = "initparms",
    nout = 1,
    outnames = c(""),
    dllname = "final"
  )

  as.data.frame(out)
}

ODEk <- function(pars) {
  x00 <- as.numeric(pars["x0"])
  y00 <- as.numeric(pars["y0"])

  rhs <- c(
    xv5 = x00, yv5 = y00,
    xv6 = x00, yv6 = y00,
    xv7 = x00, yv7 = y00,
    xv8 = x00, yv8 = y00
  )

  times <- seq(Tmin, Tmax, step_size)

  out <- ode(
    y = rhs,
    parms = pars,
    times = times,
    func = "derivs1",
    initfunc = "initparms",
    nout = 1,
    outnames = c(""),
    dllname = "final"
  )

  as.data.frame(out)
}

###################################################
## Objective function
###################################################

objective_ssr_logscale <- function(pars) {
  out <- try(ODEs(pars), silent = TRUE)
  if (inherits(out, "try-error")) return(Inf)
  if (any(!is.finite(as.matrix(out)))) return(Inf)

  ssr <- 0.0

  for (i in seq_along(fit_datatypes)) {
    ids <- obs_data[[i]][, "ids"]
    observed_values <- obs_data[[i]][, "value"]

    for (j in seq_along(ids)) {
      simulated_value <- pmax(out[ids[j], i + 1], 1e-6)
      observed_value  <- pmax(observed_values[j], 1e-6)
      ssr <- ssr + (log(simulated_value) - log(observed_value))^2
    }
  }

  ssr
}

objective_ssr_logscale_from_logpars <- function(log_pars) {
  objective_ssr_logscale(exp(log_pars))
}

sensrfun  <- function(pars) ODEs(pars)
sensrfun1 <- function(pars) ODEk(pars)

###################################################
## Parameter fitting
###################################################

lbtvec <- lbt * pars
ubtvec <- ubt * pars
lbtvec["y0"] <- 1e-3
ubtvec["y0"] <- 1e2

loglbtvec <- log(lbtvec)
logubtvec <- log(ubtvec)
cbox <- set_constr(loglbtvec, logubtvec)

theta0 <- log(pars)
theta0 <- pmin(pmax(theta0, loglbtvec + 1e-6), logubtvec - 1e-6)

est <- constrOptim(
  theta = theta0,
  f = objective_ssr_logscale_from_logpars,
  grad = NULL,
  method = "Nelder-Mead",
  ui = cbox$Amat,
  ci = cbox$bvec,
  control = list(trace = 1, maxit = 500, reltol = 1e-16)
)

fit_mixture <- exp(est$par)

print(pars)
print(objective_ssr_logscale(pars))
print(fit_mixture)
print(objective_ssr_logscale(fit_mixture))

save(fit_mixture, file = "final.Rdata")

###################################################
## Plot fitted curves
###################################################

fitted <- ODEs(pars = fit_mixture)

dd_CCC  <- fit_data_list[["condition1_ccc"]]
dd_DNA  <- fit_data_list[["condition1_dna"]]
dd_CCC1 <- fit_data_list[["condition2_ccc"]]
dd_DNA1 <- fit_data_list[["condition2_dna"]]
dd_CCC2 <- fit_data_list[["condition3_ccc"]]
dd_DNA2 <- fit_data_list[["condition3_dna"]]
dd_CCC3 <- fit_data_list[["condition4_ccc"]]
dd_DNA3 <- fit_data_list[["condition4_dna"]]

ds_CCC  <- make_sim_df(fitted, "x",   "condition1_ccc")
ds_DNA  <- make_sim_df(fitted, "y",   "condition1_dna")
ds_CCC1 <- make_sim_df(fitted, "x1s", "condition2_ccc")
ds_DNA1 <- make_sim_df(fitted, "y1s", "condition2_dna")
ds_CCC2 <- make_sim_df(fitted, "x2s", "condition3_ccc")
ds_DNA2 <- make_sim_df(fitted, "y2s", "condition3_dna")
ds_CCC3 <- make_sim_df(fitted, "x3s", "condition4_ccc")
ds_DNA3 <- make_sim_df(fitted, "y3s", "condition4_dna")

fig_plot_func(ds_CCC,  ds_DNA,  dd_CCC,  dd_DNA)
fig_plot_func(ds_CCC1, ds_DNA1, dd_CCC1, dd_DNA1)
fig_plot_func(ds_CCC2, ds_DNA2, dd_CCC2, dd_DNA2)
fig_plot_func(ds_CCC3, ds_DNA3, dd_CCC3, dd_DNA3)

###################################################
## MCMC
###################################################

# mcmc_pars <- fit_mixture
# cov0 <- 0.001 * diag(length(pars))
# MCMC_mixture <- modMCMC(
#   f = objective_ssr_logscale_from_logpars,
#   p = log(mcmc_pars),
#   lower = log(lbt * mcmc_pars),
#   upper = log(ubt * mcmc_pars),
#   niter = mcmc_iter,
#   jump = cov0,
#   var0 = NULL,
#   wvar0 = 0.1,
#   updatecov = 20
# )
# MCMC_mixture$pars <- exp(MCMC_mixture$pars)
# save(MCMC_mixture, file = "kakunin/MCMC_mixture.Rdata")

load("kakunin/MCMC_mixture.Rdata")

if (is.null(colnames(MCMC_mixture$pars))) {
  colnames(MCMC_mixture$pars) <- names(pars)[seq_len(ncol(MCMC_mixture$pars))]
}

png("kakunin/MCMC_mixture.png", width = 600, height = 600)
plot(MCMC_mixture, Full = TRUE)
pairs(MCMC_mixture, nsample = 1000)
dev.off()

cobj <- as.mcmc(MCMC_mixture$pars)
pdf("kakunin/CMC_mixture_conv.pdf")
plot(cobj)
dev.off()

###################################################
## Sensitivity ranges for fitting conditions
###################################################

sR_fit <- sensRange(func = sensrfun, parms = NULL, parInput = MCMC_mixture$pars)
save(sR_fit, file = "kakunin/MCMC_globsens.Rdata")

pdf("kakunin/MCMC_globsens.pdf")
plot(summary(sR_fit), xlab = "time")
dev.off()

###################################################
## Experimental data objects for fitted conditions
###################################################

ccc0 <- fit_data_list[["condition1_ccc"]]
dna0 <- fit_data_list[["condition1_dna"]]
ccc1 <- fit_data_list[["condition2_ccc"]]
dna1 <- fit_data_list[["condition2_dna"]]
ccc2 <- fit_data_list[["condition3_ccc"]]
dna2 <- fit_data_list[["condition3_dna"]]
ccc3 <- fit_data_list[["condition4_ccc"]]
dna3 <- fit_data_list[["condition4_dna"]]

###################################################
## Plot fitted conditions with CrI
###################################################

plot_list <- list()

fit_condition_patterns <- list(
  ccc = c("^x[[:digit:]]+\\.[[:digit:]]+$",
          "^x1s[[:digit:]]+\\.[[:digit:]]+$",
          "^x2s[[:digit:]]+\\.[[:digit:]]+$",
          "^x3s[[:digit:]]+\\.[[:digit:]]+$"),
  dna = c("^y[[:digit:]]+\\.[[:digit:]]+$",
          "^y1s[[:digit:]]+\\.[[:digit:]]+$",
          "^y2s[[:digit:]]+\\.[[:digit:]]+$",
          "^y3s[[:digit:]]+\\.[[:digit:]]+$")
)

fit_obs_pairs <- list(
  list(ccc = ccc0, dna = dna0, vline = NULL),
  list(ccc = ccc1, dna = dna1, vline = list(x = 0,  color = cmurasaki, lty = 1)),
  list(ccc = ccc2, dna = dna2, vline = list(x = 30, color = cmurasaki, lty = 1)),
  list(ccc = ccc3, dna = dna3, vline = list(x = 30, color = "black",   lty = 1))
)

for (k in seq_len(4)) {
  name_lst <- colnames(sR_fit)

  gids_ccc <- grep(fit_condition_patterns$ccc[k], name_lst)
  times_ccc <- seq(Tmin, Tmax, length = length(gids_ccc))
  c.list <- make_ci_df(sR_fit[, gids_ccc], times_ccc)

  gids_dna <- grep(fit_condition_patterns$dna[k], name_lst)
  times_dna <- seq(Tmin, Tmax, length = length(gids_dna))
  d.list <- make_ci_df(sR_fit[, gids_dna], times_dna)

  plt <- mcmc_plot_func(c.list, d.list, fit_obs_pairs[[k]]$ccc, fit_obs_pairs[[k]]$dna)

  if (!is.null(fit_obs_pairs[[k]]$vline)) {
    plt <- plt + geom_vline(
      xintercept = fit_obs_pairs[[k]]$vline$x,
      color = fit_obs_pairs[[k]]$vline$color,
      linewidth = 2,
      alpha = 0.5,
      linetype = fit_obs_pairs[[k]]$vline$lty
    )
  }

  plot_list[[k]] <- plt
  ggsave(sprintf("fig/condition%d.png", k), plt, width = 5, height = 5)
}

###################################################
## Sensitivity ranges for validation conditions
###################################################

sR_val <- sensRange(func = sensrfun1, parms = NULL, parInput = MCMC_mixture$pars)
save(sR_val, file = "kakunin/MCMC_globsens_kakunin.Rdata")

pdf("kakunin/MCMC_globsens_kakunin.pdf")
plot(summary(sR_val), xlab = "time")
dev.off()

###################################################
## Experimental data objects for validation conditions
###################################################

ccc5 <- validation_data_list[["ccc5"]]
dna5 <- validation_data_list[["dna5"]]
ccc6 <- validation_data_list[["ccc6"]]
dna6 <- validation_data_list[["dna6"]]
ccc7 <- validation_data_list[["ccc7"]]
dna7 <- validation_data_list[["dna7"]]
ccc8 <- validation_data_list[["ccc8"]]
dna8 <- validation_data_list[["dna8"]]

###################################################
## Plot validation conditions with CrI
###################################################

validation_condition_patterns <- list(
  ccc = c("^xv5[[:digit:]]+\\.[[:digit:]]+$",
          "^xv6[[:digit:]]+\\.[[:digit:]]+$",
          "^xv7[[:digit:]]+\\.[[:digit:]]+$",
          "^xv8[[:digit:]]+\\.[[:digit:]]+$"),
  dna = c("^yv5[[:digit:]]+\\.[[:digit:]]+$",
          "^yv6[[:digit:]]+\\.[[:digit:]]+$",
          "^yv7[[:digit:]]+\\.[[:digit:]]+$",
          "^yv8[[:digit:]]+\\.[[:digit:]]+$")
)

validation_obs_pairs <- list(
  list(ccc = ccc5, dna = dna5,
       vlines = list(
         list(x = 0,  color = cmurasaki, lty = 1),
         list(x = 30, color = "black",   lty = 1)
       )),
  list(ccc = ccc6, dna = dna6,
       vlines = list(
         list(x = 30, color = cmurasaki, lty = 1),
         list(x = 30, color = "black",   lty = 2)
       )),
  list(ccc = ccc7, dna = dna7,
       vlines = list(
         list(x = 30, color = cmurasaki, lty = 1),
         list(x = 45, color = "black",   lty = 1)
       )),
  list(ccc = ccc8, dna = dna8,
       vlines = list(
         list(x = 45, color = cmurasaki, lty = 1),
         list(x = 30, color = "black",   lty = 1)
       ))
)

check_C_list <- list()
check_D_list <- list()

for (k in seq_len(4)) {
  name_lst <- colnames(sR_val)

  gids_ccc <- grep(validation_condition_patterns$ccc[k], name_lst)
  times_ccc <- seq(Tmin, Tmax, length = length(gids_ccc))
  c.list <- make_ci_df(sR_val[, gids_ccc], times_ccc)

  gids_dna <- grep(validation_condition_patterns$dna[k], name_lst)
  times_dna <- seq(Tmin, Tmax, length = length(gids_dna))
  d.list <- make_ci_df(sR_val[, gids_dna], times_dna)

  plt <- mcmc_plot_func(c.list, d.list, validation_obs_pairs[[k]]$ccc, validation_obs_pairs[[k]]$dna)

  for (vl in validation_obs_pairs[[k]]$vlines) {
    plt <- plt + geom_vline(
      xintercept = vl$x,
      color = vl$color,
      linewidth = 2,
      alpha = 0.5,
      linetype = vl$lty
    )
  }

  plot_list[[k + 4]] <- plt
  ggsave(sprintf("fig/condition%d.png", k + 4), plt, width = 5, height = 5)

  check_C_list[[k]] <- make_check_df(c.list$mean, validation_obs_pairs[[k]]$ccc) %>%
    mutate(type = paste0("condition", k + 4))
  check_D_list[[k]] <- make_check_df(d.list$mean, validation_obs_pairs[[k]]$dna) %>%
    mutate(type = paste0("condition", k + 4))
}

###################################################
## Save combined figures
###################################################

fitfig <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_layout(ncol = 4)
simufig <- plot_list[[5]] + plot_list[[6]] + plot_list[[7]] + plot_list[[8]] + plot_layout(ncol = 4)

ggsave("fig/fig1.png", fitfig, width = 16, height = 4)
ggsave("fig/fig2A.png", simufig, width = 16, height = 4)

###################################################
## Fig 2B
###################################################

concat.Cdf <- bind_rows(check_C_list) %>%
  mutate(rsq = (estimated - observed)^2)

concat.Ddf <- bind_rows(check_D_list) %>%
  mutate(rsq = (estimated - observed)^2)

cdf.rsm <- sqrt(sum(concat.Cdf$rsq) / nrow(concat.Cdf))
ddf.rsm <- sqrt(sum(concat.Ddf$rsq) / nrow(concat.Ddf))

plt_ccc <- make_corr_plot(concat.Cdf, cao, "fig/fig2B1.png")
plt_dna <- make_corr_plot(concat.Ddf, caka, "fig/fig2B2.png")

###################################################
## Table 1
###################################################

print(fit_mixture)

parset <- as.data.frame(MCMC_mixture$pars)

ci_summary <- parset %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    lower  = hdi(value, credMass = 0.95)[1],
    upper  = hdi(value, credMass = 0.95)[2],
    .groups = "drop"
  )

print(ci_summary)
