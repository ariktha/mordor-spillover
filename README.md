# Geographic spillover of antimicrobial resistance (AMR) due to mass azithromycin distribution in the MORDOR trial in Niger

## Project Overview

- **Primary results of the trial:** [https://doi.org/10.1056/nejmoa1715474](https://doi.org/10.1056/nejmoa1715474) ; [https://www.nejm.org/doi/full/10.1056/NEJMc1901535](https://www.nejm.org/doi/full/10.1056/NEJMc1901535)
- **Pre-specified analysis plan:** [osf.io/u9fhc](osf.io/u9fhc)

## High-Level Script Specification

All scripts are located in `mordor-spillover/R/`

`> sessionInfo()`

`R version 4.5.0 (2025-04-11)`

`Platform: x86_64-apple-darwin20`

`Running under: macOS Sequoia 15.5`

### Core Scripts

- **00-config.R**: Project configuration, file paths, parameters, and constants.
- **00-functions.R**: Helper and core functions for distance, permutation, and analysis routines.

### Data Preparation
- **01a-exposure_dm.R**: Reads and formats treatment data from the MORDOR main (mortality) and morbidity trials.
- **01b-outcome_dm.R**: Reads and formats AMR outcome data from the MORDOR morbidity trial (24-month rectal swab DNASeq results).

### Treatment Intensity & Distance Metrics
- **02a-tx_int-idw.R**: Calculates treatment intensity at morbidity trial villages using inverse distance weighting (IDW) based on main trial treatment distribution.
- **02b-tx_int-rings.R**: Computes treatment metrics (e.g., doses) within concentric distance bands (rings) around villages.
- **02ax-layer-fns.R, 02ax-tx_int-idw-layer.R**: Additional functions and routines for spatial layer calculations.

### Correlation & Statistical Analysis
- **03a-calculate_correlations.R**: Calculates correlations between treatment intensity/doses and AMR, including permutation tests.
- **03b-plot_correlations.Rmd**: Plots permutation distributions and observed correlations.
- **03c-correlations_cis.Rmd**: Computes and visualizes confidence intervals for correlations.

### Population & Settlement Data
- **04a-crop_hrsl.R**: Crops high-resolution settlement layer (HRSL) population data to the study region.
- **04b-process_hrsl.R**: Processes HRSL data to estimate population in buffer rings around villages.
- **04c-plot_hrsl.Rmd**: Visualizes HRSL population data and buffer rings.

### Health Facility & Additional Analyses
- **05-health_posts.Rmd**: Analyzes distance to closest health facility (CSI) for each village.
- **06-log_linear_models.Rmd**: Fits log-linear models of AMR and treatment for spillover analysis.

### Manuscript Figures
- **manuscript-fig-1.Rmd, manuscript-fig-2.Rmd, manuscript-fig-3.Rmd, manuscript-supp-figs.Rmd**: Code for generating main and supplementary figures for the manuscript.

Scripts are designed to be run sequentially, starting from data preparation through to analysis and figure generation. 
**Note:** This repository does not provide geographic coordinates for trial villages or individuals, due to privacy and ethical concerns. 
