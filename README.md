# cccDNA_simulation

This repository contains the analysis code for the paper:

**“Quantification of hepatitis B virus replication dynamics in a tetracycline-off inducible cell model”**

## Overview

This repository provides the computational workflow used to analyze intracellular HBV DNA and cccDNA dynamics in a tetracycline-off inducible cell system. The code combines experimental time-course data with ordinary differential equation models to estimate model parameters, quantify uncertainty, and generate the main figures used in the manuscript.

## Repository structure

- `final.c`  
  C code defining the ODE systems used for model fitting and validation.

- `*.R`  
  R scripts for parameter estimation, MCMC-based uncertainty analysis, simulation, and figure generation.

- `data/`  
  Experimental data used for model fitting under Conditions 1–4.

- `data_kakunin/`  
  Experimental data used for validation under Conditions 5–8.

- `kakunin/`  
  Intermediate outputs such as saved MCMC objects and sensitivity-analysis results.

- `fig/`  
  Generated figures.

- `supplementary/`  

  Supplementary figures and supplementary data generated for the revised manuscript, including replicate-level measurements used in Supplementary Fig. S1 and the λ–r parameter-correlation analysis shown in Supplementary Fig. S2.
  

## Data availability

The experimental data used in this study are included in the `data/` and `data_kakunin/` folders.

## Code availability

All scripts required to reproduce the main analyses are included in this repository.

## Notes

- The model fitting for Conditions 1–4 and the validation analysis for Conditions 5–8 are handled separately in the code.
- The analysis was performed in R, with the ODE models compiled from C code for computational efficiency.
