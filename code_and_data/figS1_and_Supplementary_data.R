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
library(openxlsx)

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

        names_to = "Replicate",

        values_to = "Copies_per_well"

      ) %>%

      mutate(

        Condition = sheet,

        Marker = marker_name,

        Replicate = gsub("exp-", "Replicate ", Replicate),

        Day = as.numeric(day),

        Copies_per_well = as.numeric(Copies_per_well)

      ) %>%

      select(

        Condition,

        Marker,

        Day,

        Replicate,

        Copies_per_well

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

  filter(!is.na(Copies_per_well)) %>%

  mutate(

    Condition = factor(

      Condition,

      levels = paste0("Condition ", 1:8)

    ),

    Marker = factor(

      Marker,

      levels = c("Intracellular HBV DNA", "cccDNA")

    )

  ) %>%

  arrange(

    Condition,

    Marker,

    Day,

    Replicate

  )

###################################################

## Summary data: mean ± SD

###################################################

summary_dat <- replicate_long %>%

  group_by(

    Condition,

    Marker,

    Day

  ) %>%

  summarise(

    Mean_copies_per_well = mean(Copies_per_well, na.rm = TRUE),

    SD_copies_per_well = sd(Copies_per_well, na.rm = TRUE),

    n = sum(!is.na(Copies_per_well)),

    .groups = "drop"

  ) %>%

  arrange(

    Condition,

    Marker,

    Day

  )

###################################################

## README sheet

###################################################

readme_dat <- data.frame(

  Item = c(

    "File name",

    "Description",

    "Raw data sheet",

    "Summary sheet",

    "Condition",

    "Marker",

    "Day",

    "Replicate",

    "Copies_per_well",

    "Mean_copies_per_well",

    "SD_copies_per_well",

    "n"

  ),

  Description = c(

    "Supplementary_Data_1_replicate_level_measurements.xlsx",

    "This file contains replicate-level quantification values for intracellular HBV DNA and cccDNA used in Supplementary Fig. S1.",

    "Contains all replicate-level measurements in long format.",

    "Contains mean, standard deviation, and number of replicate wells for each condition, marker, and time point.",

    "Experimental condition corresponding to Conditions 1–8 in Fig. 1.",

    "Measured viral DNA species: Intracellular HBV DNA or cccDNA.",

    "Days post infection.",

    "Replicate well identifier.",

    "Copy number per well for each replicate measurement.",

    "Mean copy number per well across replicate wells.",

    "Standard deviation of copy number per well across replicate wells.",

    "Number of replicate wells used for summary statistics."

  ),

  stringsAsFactors = FALSE

)

###################################################

## Create Excel workbook

###################################################

supplementary_data_xlsx <- file.path(

  out_dir,

  "Supplementary_Data_1_replicate_level_measurements.xlsx"

)

wb <- createWorkbook()

addWorksheet(wb, "README")

addWorksheet(wb, "Raw data")

addWorksheet(wb, "Summary")

writeData(wb, "README", readme_dat)

writeData(wb, "Raw data", replicate_long)

writeData(wb, "Summary", summary_dat)

###################################################

## Styling

###################################################

readme_header_style <- createStyle(

  textDecoration = "bold",

  fgFill = "#E2F0D9",

  border = "Bottom"

)

header_style <- createStyle(

  textDecoration = "bold",

  fgFill = "#D9EAF7",

  border = "Bottom"

)

addStyle(

  wb,

  sheet = "README",

  style = readme_header_style,

  rows = 1,

  cols = 1:ncol(readme_dat),

  gridExpand = TRUE

)

addStyle(

  wb,

  sheet = "Raw data",

  style = header_style,

  rows = 1,

  cols = 1:ncol(replicate_long),

  gridExpand = TRUE

)

addStyle(

  wb,

  sheet = "Summary",

  style = header_style,

  rows = 1,

  cols = 1:ncol(summary_dat),

  gridExpand = TRUE

)

freezePane(wb, "README", firstRow = TRUE)

freezePane(wb, "Raw data", firstRow = TRUE)

freezePane(wb, "Summary", firstRow = TRUE)

setColWidths(wb, "README", cols = 1, widths = 28)

setColWidths(wb, "README", cols = 2, widths = 100)

setColWidths(wb, "Raw data", cols = 1:ncol(replicate_long), widths = "auto")

setColWidths(wb, "Summary", cols = 1:ncol(summary_dat), widths = "auto")

###################################################

## Save Excel file

###################################################

saveWorkbook(

  wb,

  supplementary_data_xlsx,

  overwrite = TRUE

)

cat("Saved Supplementary Data Excel file:\n")

cat(supplementary_data_xlsx, "\n")

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
  filename = file.path(out_dir, "Fig_S1_all_replicate_data_points.pdf"),
  plot = fig_s1,
  width = 18,
  height = 7
)

ggsave(
  filename = file.path(out_dir, "Fig_S1_all_replicate_data_points.png"),
  plot = fig_s1,
  width = 18,
  height = 7,
  dpi = 300
)

print(fig_s1)

###################################################
## End
###################################################
