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

library(here)
source(here("R", "00-pkg-config.R"))

run_public <- FALSE # Set to FALSE to use internal GPS data

watermark_text <- ifelse(run_public, "NOT REAL DATA", NA)
watermark_text_ml <- ifelse(run_public, "NOT\nREAL\nDATA", NA)
watermark_color <- ifelse(run_public, "darkgrey", NA)

# File paths and parallel processing ----

data_folder_path <- "/Users/ariktha/Library/CloudStorage/Box-Box/MORDOR Data/data"
path_data_clean <- here("data", "clean")
path_internal <- here("data", "internal-data")

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
recalc_idw_layer <- TRUE
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
