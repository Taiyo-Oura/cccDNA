###################################################
## Supplementary Fig. S1
## All replicate-level data points for cccDNA and intracellular HBV DNA
###################################################

rm(list = ls(all = TRUE))

###################################################
## Basic
###################################################

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(CURRENT_WORKING_DIR)

###################################################
## Packages
###################################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

###################################################
## Input files
###################################################

ccc_file <- "rawdata/cccDNA (d0-d60).xlsx"
dna_file <- "rawdata/intracellular HBV DNA (d0-d60).xlsx"

###################################################
## Output directory
###################################################

out_dir <- "supplementary"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

###################################################
## Function: read replicate-level Excel data
###################################################

read_replicate_excel <- function(file_path, marker_name) {

  sheet_names <- excel_sheets(file_path)

  dat_list <- lapply(sheet_names, function(sheet) {

    dat <- read_excel(file_path, sheet = sheet)

    dat_long <- dat %>%
      select(
        day,
        `exp-1`,
        `exp-2`,
        `exp-3`
      ) %>%
      pivot_longer(
        cols = starts_with("exp-"),
        names_to = "replicate",
        values_to = "copies_per_well"
      ) %>%
      mutate(
        condition = sheet,
        marker = marker_name,
        replicate = gsub("exp-", "Replicate ", replicate),
        day = as.numeric(day),
        copies_per_well = as.numeric(copies_per_well)
      )

    dat_long
  })

  bind_rows(dat_list)
}

###################################################
## Read data
###################################################

ccc_long <- read_replicate_excel(
  file_path = ccc_file,
  marker_name = "cccDNA"
)

dna_long <- read_replicate_excel(
  file_path = dna_file,
  marker_name = "Intracellular HBV DNA"
)

replicate_long <- bind_rows(ccc_long, dna_long) %>%
  filter(!is.na(copies_per_well)) %>%
  mutate(
    condition = factor(
      condition,
      levels = paste0("Condition ", 1:8)
    ),
    marker = factor(
      marker,
      levels = c("Intracellular HBV DNA", "cccDNA")
    )
  )

###################################################
## Summary data: mean ± SD
## Saved only as supplementary data, not shown in Fig. S1
###################################################

summary_dat <- replicate_long %>%
  group_by(condition, marker, day) %>%
  summarise(
    mean = mean(copies_per_well, na.rm = TRUE),
    sd = sd(copies_per_well, na.rm = TRUE),
    n = sum(!is.na(copies_per_well)),
    .groups = "drop"
  )

###################################################
## Save Supplementary Data CSV
###################################################

write.csv(
  replicate_long,
  file = file.path(out_dir, "Supplementary_Data_1_replicate_level_measurements_long.csv"),
  row.names = FALSE
)

write.csv(
  summary_dat,
  file = file.path(out_dir, "Supplementary_Data_1_replicate_level_measurements_summary.csv"),
  row.names = FALSE
)

###################################################
## Plot function
###################################################

plot_marker <- function(marker_name) {

  plot_dat <- replicate_long %>%
    filter(marker == marker_name)

  ggplot(plot_dat, aes(x = day, y = copies_per_well)) +
    geom_jitter(
      width = 0,
      height = 0,
      size = 1.2,
      stroke = 0.4,
      alpha = 1,
      shape = 21,
      color = "black",
      fill = NA
    ) +
    facet_wrap(~ condition, nrow = 1) +
    scale_y_log10(
      limits = c(1e0, 1e12),
      breaks = c(1e0, 1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
      labels = c(
        expression(10^0),
        expression(10^2),
        expression(10^4),
        expression(10^6),
        expression(10^8),
        expression(10^10),
        expression(10^12)
      )
    ) +
    scale_x_continuous(
      limits = c(0, 60),
      breaks = seq(0, 60, by = 20)
    ) +
    labs(
      title = marker_name,
      x = "Days post infection",
      y = "Copies/well"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      legend.position = "none"
    )
}

###################################################
## Make Supplementary Fig. S1
###################################################

p_dna <- plot_marker("Intracellular HBV DNA")
p_ccc <- plot_marker("cccDNA")

fig_s1 <- p_dna / p_ccc

###################################################
## Save figure
###################################################

ggsave(
  filename = file.path(out_dir, "Supplementary_Fig_S1_all_replicate_data_points.pdf"),
  plot = fig_s1,
  width = 18,
  height = 7
)

ggsave(
  filename = file.path(out_dir, "Supplementary_Fig_S1_all_replicate_data_points.png"),
  plot = fig_s1,
  width = 18,
  height = 7,
  dpi = 300
)

print(fig_s1)

###################################################
## End
###################################################
