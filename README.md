# cccDNA_simulation

This repository contains the analysis code and data associated with the manuscript:

**“Quantification of hepatitis B virus replication dynamics in a tetracycline-off inducible cell model”**

## Overview

This repository provides the computational workflow used to analyze intracellular HBV DNA and covalently closed circular DNA (cccDNA) dynamics in a tetracycline-off inducible HBV cell system.

The workflow combines experimental time-course measurements with ordinary differential equation models to:

- estimate model parameters using data from Conditions 1–4,
- assess parameter uncertainty using Markov chain Monte Carlo sampling,
- validate the model using independent data from Conditions 5–8,
- calculate uncertainty and predictive intervals,
- examine parameter relationships and practical identifiability, and
- generate the figures included in the manuscript and Supplementary Information.

All analysis code and data are located in the `code_and_data/` directory.

## Repository structure

- `code_and_data/final.c`  
  C source code defining the ordinary differential equation systems used for model fitting and validation.

- `code_and_data/*.R`  
  R scripts for parameter estimation, MCMC-based uncertainty analysis, model simulation, validation, sensitivity analysis, parameter-correlation analysis, and figure generation.

- `code_and_data/data/`  
  Experimental time-course data used for model fitting under Conditions 1–4.

- `code_and_data/data_kakunin/`  
  Independent experimental time-course data used for model validation under Conditions 5–8.

- `code_and_data/kakunin/`  
  Intermediate analysis outputs, including saved MCMC objects and sensitivity-analysis results.

- `code_and_data/fig/`  
  Figures generated from the main analysis scripts.

- `code_and_data/rawdata/`  
  Replicate-level measurements of intracellular HBV DNA and cccDNA.

- `code_and_data/supplementary/`  
  Scripts, figures, and data associated with the Supplementary Information. This folder includes:

  - replicate-level measurements and plots used in Fig. S1, and
  - predictive-interval analyses used in Fig. S2.

## Analysis workflow

### 1. Model fitting

Conditions 1–4 were used for simultaneous model fitting.

Model parameters were estimated by minimizing the residual sum of squares between the log-transformed model predictions and experimental measurements.

### 2. Parameter uncertainty

Parameter uncertainty was assessed using adaptive Metropolis–Hastings Markov chain Monte Carlo sampling.

The MCMC parameter samples were propagated through the ordinary differential equation model to generate trajectory distributions and 95% MCMC-derived uncertainty intervals.

### 3. Predictive intervals

Residual observation error was estimated from the log-scale residuals between the experimental measurements and the best-fit model predictions for Conditions 1–4.

This residual error was added to trajectories generated from the MCMC parameter samples to construct the 95% predictive intervals shown in Fig. S2.

### 4. Model validation

Conditions 5–8 were used as independent validation datasets.

No additional parameter fitting was performed for these conditions. The parameter estimates obtained from Conditions 1–4 were used to simulate the corresponding treatment schedules.

### 5. Parameter-correlation analysis

The relationship between the empirical parameters λ and r was evaluated using MCMC-sampled parameter values and a two-dimensional profile residual-sum-of-squares analysis.

These analyses were used to assess their practical identifiability and compensatory relationship.

## Running the analysis

The analysis scripts are intended to be run in RStudio.

Open the repository in RStudio, set `code_and_data/` as the working directory, and run the corresponding R scripts.

For example:

```r
setwd("code_and_data")
source("figS2.R")

```

Some scripts determine their working directory automatically using:

```r
CURRENT_WORKING_DIR <- dirname(
  rstudioapi::getActiveDocumentContext()$path
)
setwd(CURRENT_WORKING_DIR)
```

Therefore, these scripts should be opened and executed directly in RStudio.

The C source file `final.c` is compiled from R within the relevant analysis scripts. Compiled object and shared-library files are not included in the repository because they are generated locally and may depend on the operating system and R environment.

## Data availability

The experimental time-course data used for model fitting and validation are included in the `code_and_data/data/` and `code_and_data/data_kakunin/` folders.

Replicate-level raw measurements are provided in the `code_and_data/rawdata/` folder.

Supplementary figures and associated processed data are provided in the `code_and_data/supplementary/` folder.

## Code availability

All scripts required to reproduce the main analyses, validation analyses, uncertainty calculations, predictive intervals, and associated figures are included in this repository.

## Requirements

The analysis was performed in R.

The ordinary differential equation models were implemented in C and compiled from R for computational efficiency.

The main R packages used in the analysis include:

- `deSolve`
- `FME`
- `ggplot2`
- `dplyr`
- `tidyr`
- `purrr`
- `patchwork`
- `rstudioapi`

A working C compiler is required to compile `final.c`.

## Notes

- Conditions 1–4 were used for model fitting.
- Conditions 5–8 were used for independent validation.
- No parameter refitting was performed for Conditions 5–8.
- Model fitting was performed on the log-transformed scale.
- The MCMC-derived uncertainty intervals represent uncertainty arising from parameter estimation.
- The predictive intervals additionally incorporate residual observation error estimated from Conditions 1–4.
- The ODE models are compiled from `final.c` before model simulation.

## License

The license for this repository is provided in the `LICENSE` file.
