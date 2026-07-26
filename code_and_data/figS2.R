###################################################

rm(list = ls(all = TRUE))

###################################################
## Basic
###################################################

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(CURRENT_WORKING_DIR)

###################################################
## Package
###################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(deSolve)

###################################################
## Load
###################################################

load("final.Rdata")
load("kakunin/MCMC_mixture.Rdata")

###################################################
## Output directory
###################################################

dir.create("supplementary", showWarnings = FALSE)

###################################################
## Compile ODE model
###################################################

system("R CMD SHLIB final.c")
dyn.load(paste0("final", .Platform$dynlib.ext))

###################################################
## Settings
###################################################

Tmin <- 0
Tmax <- 60
step_size <- 0.1
stime <- seq(Tmin, Tmax, step_size)

set.seed(1234)

N_MCMC_DRAW <- 2000
MIN_VALUE <- 1e-12

###################################################
## Experimental data
###################################################

read_txt_data <- function(file_name) {
  dat <- read.csv(
    file_name,
    sep = "\t",
    header = FALSE
  )

  data.frame(
    time = dat[, 1],
    value = dat[, 2]
  )
}

read_csv_data <- function(file_name) {
  dat <- read.csv(
    file_name,
    header = FALSE
  )

  data.frame(
    time = dat[, 1],
    value = dat[, 2]
  )
}

observed_data <- bind_rows(

  lapply(
    1:4,
    function(i) {
      bind_rows(

        read_txt_data(
          paste0("data/condition", i, "_ccc.txt")
        ) %>%
          mutate(
            condition = i,
            variable = "cccDNA"
          ),

        read_txt_data(
          paste0("data/condition", i, "_dna.txt")
        ) %>%
          mutate(
            condition = i,
            variable = "Intracellular HBV DNA"
          )
      )
    }
  ),

  lapply(
    5:8,
    function(i) {
      bind_rows(

        read_csv_data(
          paste0("data_kakunin/ccc", i, ".csv")
        ) %>%
          mutate(
            condition = i,
            variable = "cccDNA"
          ),

        read_csv_data(
          paste0("data_kakunin/dna", i, ".csv")
        ) %>%
          mutate(
            condition = i,
            variable = "Intracellular HBV DNA"
          )
      )
    }
  )
) %>%
  mutate(
    value = pmax(value, MIN_VALUE)
  )

###################################################
## ODE functions
###################################################

ODEs <- function(pars) {

  x00 <- pars["x0"]
  y00 <- pars["y0"]

  initial_values <- c(
    x = x00,
    y = y00,
    x1s = x00,
    y1s = y00,
    x2s = x00,
    y2s = y00,
    x3s = x00,
    y3s = y00
  )

  out <- ode(
    y = initial_values,
    parms = pars,
    times = stime,
    func = "derivs",
    initfunc = "initparms",
    nout = 1,
    outnames = c(""),
    dllname = "final"
  )

  as.data.frame(out)
}

ODEk <- function(pars) {

  x00 <- pars["x0"]
  y00 <- pars["y0"]

  initial_values <- c(
    xv5 = x00,
    yv5 = y00,
    xv6 = x00,
    yv6 = y00,
    xv7 = x00,
    yv7 = y00,
    xv8 = x00,
    yv8 = y00
  )

  out <- ode(
    y = initial_values,
    parms = pars,
    times = stime,
    func = "derivs1",
    initfunc = "initparms",
    nout = 1,
    outnames = c(""),
    dllname = "final"
  )

  as.data.frame(out)
}

simulate_all_conditions <- function(pars) {

  fit_out <- ODEs(pars)
  val_out <- ODEk(pars)

  bind_rows(

    data.frame(
      condition = 1,
      time = fit_out$time,
      cccDNA = fit_out$x,
      HBV_DNA = fit_out$y
    ),

    data.frame(
      condition = 2,
      time = fit_out$time,
      cccDNA = fit_out$x1s,
      HBV_DNA = fit_out$y1s
    ),

    data.frame(
      condition = 3,
      time = fit_out$time,
      cccDNA = fit_out$x2s,
      HBV_DNA = fit_out$y2s
    ),

    data.frame(
      condition = 4,
      time = fit_out$time,
      cccDNA = fit_out$x3s,
      HBV_DNA = fit_out$y3s
    ),

    data.frame(
      condition = 5,
      time = val_out$time,
      cccDNA = val_out$xv5,
      HBV_DNA = val_out$yv5
    ),

    data.frame(
      condition = 6,
      time = val_out$time,
      cccDNA = val_out$xv6,
      HBV_DNA = val_out$yv6
    ),

    data.frame(
      condition = 7,
      time = val_out$time,
      cccDNA = val_out$xv7,
      HBV_DNA = val_out$yv7
    ),

    data.frame(
      condition = 8,
      time = val_out$time,
      cccDNA = val_out$xv8,
      HBV_DNA = val_out$yv8
    )
  ) %>%
    mutate(
      cccDNA = pmax(cccDNA, MIN_VALUE),
      HBV_DNA = pmax(HBV_DNA, MIN_VALUE)
    )
}

###################################################
## Best-fit trajectories
###################################################

best_fit_trajectory <- simulate_all_conditions(
  fit_mixture
)

###################################################
## Residual observation error
## Conditions 1–4 only
###################################################

best_fit_long <- best_fit_trajectory %>%
  pivot_longer(
    cols = c(
      cccDNA,
      HBV_DNA
    ),
    names_to = "variable_code",
    values_to = "predicted"
  ) %>%
  mutate(
    variable = ifelse(
      variable_code == "cccDNA",
      "cccDNA",
      "Intracellular HBV DNA"
    )
  )

residual_data <- observed_data %>%
  filter(
    condition %in% 1:4
  ) %>%
  group_by(
    condition,
    variable
  ) %>%
  group_modify(
    function(.x, .y) {

      model_data <- best_fit_long %>%
        filter(
          condition == .y$condition,
          variable == .y$variable
        )

      .x %>%
        mutate(
          predicted = approx(
            x = model_data$time,
            y = model_data$predicted,
            xout = time,
            rule = 2
          )$y,

          residual =
            log(value) -
            log(predicted)
        )
    }
  ) %>%
  ungroup()

sigma_cccDNA <- residual_data %>%
  filter(
    variable == "cccDNA"
  ) %>%
  summarise(
    sigma = sd(residual)
  ) %>%
  pull(sigma)

sigma_HBV_DNA <- residual_data %>%
  filter(
    variable == "Intracellular HBV DNA"
  ) %>%
  summarise(
    sigma = sd(residual)
  ) %>%
  pull(sigma)

###################################################
## MCMC samples
###################################################

parameter_samples <- as.data.frame(
  MCMC_mixture$pars
)

parameter_samples <- parameter_samples %>%
  filter(
    if_all(
      everything(),
      ~ is.finite(.x) & .x > 0
    )
  )

if (nrow(parameter_samples) > N_MCMC_DRAW) {
  parameter_samples <- parameter_samples %>%
    slice_sample(
      n = N_MCMC_DRAW
    )
}

###################################################
## MCMC trajectories
###################################################

mcmc_trajectories <- map_dfr(
  1:nrow(parameter_samples),
  function(i) {

    pars_i <- unlist(
      parameter_samples[i, ]
    )

    simulate_all_conditions(
      pars_i
    ) %>%
      mutate(
        draw_id = i
      )
  }
)

###################################################
## Predictive intervals
###################################################

predictive_samples <- mcmc_trajectories %>%
  mutate(
    cccDNA_predictive = exp(
      log(cccDNA) +
        rnorm(
          n(),
          mean = 0,
          sd = sigma_cccDNA
        )
    ),

    HBV_DNA_predictive = exp(
      log(HBV_DNA) +
        rnorm(
          n(),
          mean = 0,
          sd = sigma_HBV_DNA
        )
    )
  )

interval_data <- predictive_samples %>%
  group_by(
    condition,
    time
  ) %>%
  summarise(
    cccDNA_median =
      median(
        cccDNA,
        na.rm = TRUE
      ),

    cccDNA_lower =
      quantile(
        cccDNA_predictive,
        0.025,
        na.rm = TRUE
      ),

    cccDNA_upper =
      quantile(
        cccDNA_predictive,
        0.975,
        na.rm = TRUE
      ),

    HBV_DNA_median =
      median(
        HBV_DNA,
        na.rm = TRUE
      ),

    HBV_DNA_lower =
      quantile(
        HBV_DNA_predictive,
        0.025,
        na.rm = TRUE
      ),

    HBV_DNA_upper =
      quantile(
        HBV_DNA_predictive,
        0.975,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

###################################################
## Plot data
###################################################

plot_data <- bind_rows(

  interval_data %>%
    transmute(
      condition,
      time,
      variable = "cccDNA",
      median = cccDNA_median,
      lower = cccDNA_lower,
      upper = cccDNA_upper
    ),

  interval_data %>%
    transmute(
      condition,
      time,
      variable = "Intracellular HBV DNA",
      median = HBV_DNA_median,
      lower = HBV_DNA_lower,
      upper = HBV_DNA_upper
    )
) %>%
  mutate(
    condition_label = factor(
      paste0(
        "Condition ",
        condition
      ),
      levels = paste0(
        "Condition ",
        1:8
      )
    )
  )

observed_plot_data <- observed_data %>%
  mutate(
    condition_label = factor(
      paste0(
        "Condition ",
        condition
      ),
      levels = paste0(
        "Condition ",
        1:8
      )
    )
  )

###################################################
## Fig. S2
###################################################

Fig_S2 <- ggplot() +

  geom_ribbon(
    data = plot_data,
    aes(
      x = time,
      ymin = lower,
      ymax = upper,
      fill = variable
    ),
    alpha = 0.20,
    color = NA
  ) +

  geom_line(
    data = plot_data,
    aes(
      x = time,
      y = median,
      color = variable
    ),
    linewidth = 0.9
  ) +

  geom_point(
    data = observed_plot_data,
    aes(
      x = time,
      y = value,
      color = variable
    ),
    size = 1.6
  ) +

  facet_wrap(
    ~ condition_label,
    ncol = 4
  ) +

  scale_y_log10(
    limits = c(
      1,
      1e12
    )
  ) +

  scale_x_continuous(
    limits = c(
      0,
      60
    ),
    breaks = seq(
      0,
      60,
      10
    )
  ) +

  scale_color_manual(
    values = c(
      "cccDNA" = "#3366CC",
      "Intracellular HBV DNA" = "#CC3366"
    )
  ) +

  scale_fill_manual(
    values = c(
      "cccDNA" = "#3366CC",
      "Intracellular HBV DNA" = "#CC3366"
    )
  ) +

  labs(
    x = "Days post infection",
    y = "cccDNA / intracellular HBV DNA (copies/well)",
    color = NULL,
    fill = NULL
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.position = "bottom"
  )

print(Fig_S2)

###################################################
## Save
###################################################

ggsave(
  "supplementary/Fig_S2_predictive_intervals.png",
  Fig_S2,
  width = 12,
  height = 6,
  dpi = 600
)

ggsave(
  "supplementary/Fig_S2_predictive_intervals.pdf",
  Fig_S2,
  width = 12,
  height = 6
)
