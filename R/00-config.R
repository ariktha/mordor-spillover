#--------------------------------------------
#
# mordor-spillover/R/00-config.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Config file
#
#------------------------------------

rm(list = ls())

library(tidyverse)

library(sf)
sf_use_s2(FALSE)

library(parallel)

data_folder_path <- "/Users/ariktha/Library/CloudStorage/Box-Box/MORDOR Data/data"
n_cores <- max(1, parallel::detectCores() - 2)

# Coordinate Reference Systems ----

## Primarily use Global CRS: WGS 84
global_crs <- 4326

## Projected CRS: UTM zone 31N
## Coordinates are projected to UTM zone 31N (specific to Niger) for distance calculations
proj_crs <- 32631

# Parameters for exposure measures ----

## IDW parameters
idw_power <- 1
trunc_val <- c(1, 1)

## IDW layer re-estimation
recalc_idw_layer <- FALSE
pred_grid_dim <- 300

## Radii for doses within distance
primary_out_radii <- seq(10, 50, by = 10)

# Antimicrobial resistance classes of interest ----

ab_classes_of_interest <- c("Aminocoumarins", "Aminoglycosides", "Bacitracin", "betalactams", "Elfamycins", 
                            "Fluoroquinolones", "Fosfomycin", "Glycopeptides", "Metronidazole", "MLS",
                            "Multi.drug.resistance", "Phenicol", "Rifampin", "Sulfonamides",
                            "Sulfonamides", "Tetracyclines", "Trimethoprim")

# Permutation test ----

# Number of permutations for significance testing
n_permutations <- 1000

# Color palette for AMR classes ----

arms_colors <- c("azithro" = "#f1a226", "placebo" = "#003f5c")
