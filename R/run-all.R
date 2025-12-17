renv::restore()

library(here)
library(quarto)
library(knitr)

# Download public data
source(here("R", "00-download-public-data.R"))

if (!dir.exists(path_data_clean)) {
  dir.create(here("R", "markdown"), recursive = TRUE)
}

# Primary analysis: calculate exposure measures and run models
source(here("R", "02a-tx_int-idw.R"))
source(here("R", "02b-tx_int-rings.R"))

source(here("R", "03-calculate_correlations.R"))

# Potential confounders: 

## High Resolution Settlement Layer (HRSL) population density
source(here("R", "04a-crop_hrsl.R"))
source(here("R", "04b-process_hrsl.R"))

## Distance to nearest CSI
source(here("R", "05a-process_csi_dist.R"))
knit(here("R", "05b-morb_v_main_csi.Rmd"), 
     output = here("R", "markdown", "05b-morb_v_main_csi.md"))

# Log-linear models and change in resistance analysis
knit(here("R", "06-log_linear_models.Rmd"), 
     output = here("R", "markdown", "06-log_linear_models.md"))
knit(here("R", "07-change-in-mls-analysis.Rmd"), 
     output = here("R", "markdown", "07-change-in-mls-analysis.md"))

# Individual-level analysis
quarto_render(here("R", "08-indiv-data.qmd"))

# Phenotypic resistance analysis
quarto_render(here("R", "09-mordor-spillover-NP.qmd"))

# Figures for manuscript
knit(here("R", "manuscript-fig-1.Rmd"), 
     output = here("R", "markdown", "manuscript-fig-1.md"))
knit(here("R", "manuscript-fig-2.Rmd"), 
     output = here("R", "markdown", "manuscript-fig-2.md"))
knit(here("R", "manuscript-fig-3.Rmd"), 
     output = here("R", "markdown", "manuscript-fig-3.md"))

## Supplementary figures for manuscript
knit(here("R", "manuscript-supp-figs.Rmd"), 
     output = here("R", "markdown", "manuscript-supp-figs.md"))