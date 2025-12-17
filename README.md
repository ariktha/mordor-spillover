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

If you have any questions about the files in this repository, please contact Ariktha Srivathsan or Ben Arnold at UCSF (ariktha.srivathsan@ucsf.edu; ben.arnold@ucsf.edu)

## Code Submission Items

Following: https://www.nature.com/documents/nr-software-policy.pdf

### System requirements

`> sessionInfo()`

`R version 4.5.1 (2025-06-13)`

`Platform: x86_64-apple-darwin20`

`Running under: macOS Sequoia 15.6.1`

All analyses were run using R version 4.5.1 (2025-06-13) -- "Great Square Root" on Mac OSX Sequoia (15.6.1) using the RStudio IDE (https://www.rstudio.com). Analyses were also tested on R v4.4.0.
In this repository we have used the `renv` package to archive the package versions so that you and reproduce the exact compute environment, should you wish to do so. 

### Installation Guide and Instructions for Use

- You can download and install R from CRAN: https://cran.r-project.org
- You can download and install RStudio from their website: https://www.rstudio.com
- All R packages required to run the analyses are sourced in the file `00-pkg-config.R`.
- The installation time should be < 10 minutes total on a typical desktop computer.

To reproduce all analyses in the paper, we recommend that you: 

1. Clone the GitHub repository to your computer using `git clone https://github.com/ariktha/mordor-spillover.git` in the terminal

2. Recreate the exact package environment using the `renv` package. 

 You can do this by opening the R project file ([mordor-spillover.Rproj](https://github.com/ariktha/mordor-spillover/blob/main/mordor-spillover.Rproj)) in RStudio, loading the `renv` package, and typing `renv::restore()` to restore the package environment from the projects [renv.lock](https://github.com/ariktha/mordor-spillover/blob/main/renv.lock) file. 

3. Download the public data from the OSF repository by running the script [`00-download-public-data.R`](https://github.com/ariktha/mordor-spillover/blob/main/R/00-download-public-data.R).
  
4. All of the analysis scripts can be run sequentially using the script [`run-all.R`](https://github.com/ariktha/mordor-spillover/blob/main/R/run-all.R).

### Additional details

- Scripts `01a-exposure_dm.R` and `01b-outcome_dm.R` process raw data and are provided for completeness. The resulting datasets are shared publicly in https://osf.io/bmjd3/ in csv and rds formats along with their codebooks.
- These datasets can be downloaded using the script `00-download-public-data.R`.
- GPS coordinates for village and health post locations cannot be made public due to privacy considerations. To facilitate code reproducibility, `00-download-public-data.R` generates synthetic example datasets with the same structure as the original data. Results produced using these synthetic datasets will differ from published results, which are based on the actual data.

## Repository Structure

### Core Scripts

- **00-config.R**: Configuration file with global parameters, file paths, and study constants.
- **00-pkg-config.R**: Loads all required packages.
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
- **05c-morb_v_main_csi.Rmd**

### Modelling
- **06-log_linear_models.Rmd**: Fits log-linear models of AMR and treatment for spillover analysis.
- **07-change-in-mls-analysis.Rmd**: Estimates correlations between the change in AMR reads and geographic treatment intensity.

### Supplementary analyses
- **08-indiv-data.qmd**
- **09-mordor-spillover-NP.qmd**

### Manuscript Figures
- **manuscript-fig-1.Rmd, manuscript-fig-2.Rmd, manuscript-fig-3.Rmd, manuscript-supp-figs.Rmd**: Code for generating main and supplementary figures for the manuscript.

## Data Privacy Note
GPS coordinates for village locations cannot be made public due to privacy considerations. Additionally, since each permutation test re-estimates the treatment intensity layer using randomized treatment assignments, the complete spatial mapping scripts and statistical inference (p-values) cannot be reproduced using the public data alone.
