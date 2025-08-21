# mordor-spillover

## Geographic spillover of antimicrobial resistance (AMR) due to mass azithromycin distribution in the MORDOR trial in Niger

## Project Overview

Large-scale, placebo-controlled, cluster-randomized trials in high-mortality settings in several African countries demonstrated a 14-18% reduction in childhood mortality following twice-annual mass drug administration (MDA) of azithromycin among children aged 1–59 months. Azithromycin MDA also selects for antimicrobial resistance (AMR), particularly macrolide resistance in treated populations. It is unknown whether the genetic selection of AMR from azithromycin MDA could spill over to neighboring untreated populations.
We assessed between-village geographic spillover effects of genotypic resistance to macrolides and other antibiotic classes in rectal swabs collected from 1200 children in 30 monitoring villages in Niger after two years of MDA in 594 surrounding villages.
Key Finding: We found no evidence of geographic spillover of macrolide resistance in untreated villages, as the genetic load of AMR remained at baseline levels in placebo-treated villages regardless of surrounding azithromycin treatment intensity.

- **Preprint:** [www.medrxiv.org/content/10.1101/2025.07.22.25331994v1](https://www.medrxiv.org/content/10.1101/2025.07.22.25331994v1)
- **Pre-specified analysis plan:** [osf.io/u9fhc](osf.io/u9fhc)
- **Public data repository:** [osf.io/bmjd3](https://osf.io/bmjd3/?view_only=209cd7c09cb44fa3bd76e9a61e5f805e)
- **Primary results of the trials:**
  - MORDOR Mortality trial: [https://doi.org/10.1056/nejmoa1715474](https://doi.org/10.1056/nejmoa1715474)
  - MORDOR Morbidity (AMR-monitoring) trial: [https://www.nejm.org/doi/full/10.1056/NEJMc1901535](https://www.nejm.org/doi/full/10.1056/NEJMc1901535)


## High-Level Script Specification

All scripts are located in `mordor-spillover/R/`

`> sessionInfo()`

`R version 4.5.0 (2025-04-11)`

`Platform: x86_64-apple-darwin20`

`Running under: macOS Sequoia 15.5`

## Repository Structure

### Core Scripts

- **00-config.R**: Configuration file with global parameters, file paths, and study constants.
- **00-functions.R**: Helper and custom functions for distance calculations, permutation tests, and analysis routines.

### Data Preparation
- **01a-exposure_dm.R**: Reads and formats treatment data from the MORDOR mortality and morbidity (AMR-monitoring) trials.
- **01b-outcome_dm.R**: Reads and formats AMR outcome data from the MORDOR morbidity trial (24-month rectal swab DNASeq results).

### Treatment Intensity & Distance Metrics
- **02a-tx_int-idw.R**: Calculates geographic treatment intensity at morbidity trial villages using inverse distance weighting based on mortality trial treatment distribution.
- **02b-tx_int-rings.R**: Computes cumulative treatment metrics within concentric distance bands (rings) around AMR-monitoring villages.
- **02ax-layer-fns.R, 02ax-tx_int-idw-layer.R**: Functions and routine to generate a geographic treatment intensity raster layer for visualization.

### Correlation & Statistical Analysis
- **03a-calculate_correlations.R**: Calculates correlations between treatment intensity/doses and AMR with permutation tests for inference.
- **03b-plot_correlations.Rmd**: Plots permutation distributions and observed correlations.
- **03c-correlations_cis.Rmd**: Computes and visualizes confidence intervals for correlations.

### Confounding assessment: Population & Settlement Data
- **04a-crop_hrsl.R**: Crops high-resolution settlement layer (HRSL) population data to the study region.
- **04b-process_hrsl.R**: Processes HRSL data to estimate population in buffer rings around villages.
- **04c-plot_hrsl.Rmd**: Visualizes HRSL population data and buffer rings.

### Confounding assessment: Health Facility data
- **05-health_posts.Rmd**: Analyzes distance to closest health facility (CSI) for each village.

### Modelling
- **06-log_linear_models.Rmd**: Fits log-linear models of AMR and treatment for spillover analysis.

### Manuscript Figures
- **manuscript-fig-1.Rmd, manuscript-fig-2.Rmd, manuscript-fig-3.Rmd, manuscript-supp-figs.Rmd**: Code for generating main and supplementary figures for the manuscript.

## Data Privacy Note
GPS coordinates for village locations cannot be made public due to privacy considerations. Additionally, since each permutation test re-estimates the treatment intensity layer using randomized treatment assignments, the complete spatial mapping scripts and statistical inference (p-values) cannot be reproduced using the public data alone.
