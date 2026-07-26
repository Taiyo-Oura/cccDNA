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

## Repository structure

- `final.c`  
  C source code defining the ordinary differential equation systems used for model fitting and validation.

- `*.R`  
  R scripts for parameter estimation, MCMC-based uncertainty analysis, model simulation, validation, sensitivity analysis, parameter-correlation analysis, and figure generation.

- `data/`  
  Experimental time-course data used for model fitting under Conditions 1–4.

- `data_kakunin/`  
  Independent experimental time-course data used for model validation under Conditions 5–8.

- `kakunin/`  
  Intermediate analysis outputs, including saved MCMC objects and sensitivity-analysis results.

- `fig/`  
  Figures generated from the main analysis scripts.

- `rawdata/`  
  Replicate-level measurements of intracellular HBV DNA and cccDNA.

- `supplementary/`  
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

## Data availability

The experimental time-course data used for model fitting and validation are included in the `data/` and `data_kakunin/` folders.

Replicate-level raw measurements are provided in the `rawdata/` folder.

Supplementary figures and associated processed data are provided in the `supplementary/` folder.

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

## Notes

- Conditions 1–4 were used for model fitting.
- Conditions 5–8 were used for independent validation.
- No parameter refitting was performed for Conditions 5–8.
- Model fitting was performed on the log-transformed scale.
- The MCMC-derived uncertainty intervals represent uncertainty arising from parameter estimation.
- The predictive intervals additionally incorporate residual observation error estimated from Conditions 1–4.
- The ODE models are compiled from `final.c` before model simulation.
