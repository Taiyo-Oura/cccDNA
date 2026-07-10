# cccDNA_simulation

This repository contains the analysis code and data associated with the manuscript:

**“Quantification of hepatitis B virus replication dynamics in a tetracycline-off inducible cell model”**

## Overview

This repository provides the computational workflow used to analyze intracellular HBV DNA and cccDNA dynamics in a tetracycline-off inducible cell system. The workflow combines experimental time-course data with ordinary differential equation models to estimate model parameters, quantify parameter uncertainty, perform validation analyses, and generate the figures used in the manuscript and Supplementary Information.

## Repository structure

- `final.c`  
  C source code defining the ordinary differential equation systems used for model fitting and validation.

- `*.R`  
  R scripts for parameter estimation, MCMC-based uncertainty analysis, simulation, sensitivity analysis, and figure generation.

- `data/`  
  Experimental time-course data used for model fitting under Conditions 1–4.

- `data_kakunin/`  
  Experimental time-course data used for validation under Conditions 5–8.

- `kakunin/`  
  Intermediate outputs, including saved MCMC objects and sensitivity-analysis results.

- `fig/`  
  Figures generated from the main analysis scripts.

- `rawdata/`  
  Replicate-level raw measurement data for intracellular HBV DNA and cccDNA.

- `supplementary/`  
  Supplementary figures and supplementary data generated for the revised manuscript. This folder includes the replicate-level measurements used in Supplementary Fig. S1 and the λ–r parameter-correlation analysis shown in Supplementary Fig. S2.

## Data availability

The experimental time-course data used for model fitting and validation are included in the `data/` and `data_kakunin/` folders. Replicate-level raw measurement data are provided in the `rawdata/` folder, and processed supplementary data are provided in the `supplementary/` folder.

## Code availability

All scripts required to reproduce the main analyses and generate the associated figures are included in this repository.

## Notes

- Conditions 1–4 were used for model fitting.
- Conditions 5–8 were used for validation.
- The analysis was performed in R.
- The ODE models were implemented in C and compiled from R for computational efficiency.
