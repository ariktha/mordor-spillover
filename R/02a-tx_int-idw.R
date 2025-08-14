#--------------------------------------------
#
# mordor-spillover/R/02a-tx_int-idw.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Calculate treatment intensity at morbidity trial villages based on
# the main trial treatment distribution.
#
#------------------------------------

library(here)

rm(list = ls())       

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))

# Load data ---------------------------------------------------------------

# Load MORDOR morbidity and main trial data

morb_grappe <- readRDS(here("data", "clean", "morb_grappe.RDS"))
main_grappe <- readRDS(here("data", "clean", "main_grappe.RDS"))

# Load gps data

morb_gps <- readRDS(here("data", "clean", "morb_gps.RDS"))
main_gps <- readRDS(here("data", "clean", "main_gps.RDS"))

# Create distance matrix --------------------------------------------------

# Prepare distance matrix for morbidity and main trial data

morb_grappe <- morb_grappe %>% 
  left_join(morb_gps, by = "grappe") %>%
  rename(morbidity_grappe = grappe) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = global_crs, remove = FALSE)

main_grappe <- main_grappe %>% 
  dplyr::select(grappe, arm, treat_bin, n_doses, n_children) %>%
  left_join(main_gps, by = "grappe") %>%
  rename(main_grappe = grappe) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = global_crs, remove = FALSE)

dist_main_morb <- compute_distance_matrix(
  df1 = main_grappe, 
  df2 = morb_grappe, 
  colname1 = "main_grappe", 
  colname2 = "morbidity_grappe") %>%
  left_join(main_grappe %>% st_drop_geometry(), 
            by = "main_grappe")

saveRDS(dist_main_morb, here("data", "output", "dist_main_morb.RDS"))

# Calculate treatment intensity at morbidity trial villages ---------------

morb_tx_int <- calculate_idw(
  dist_long = dist_main_morb, 
  measure_at_var = "morbidity_grappe", 
  numerator_var = "n_doses",
  power = c(idw_power, 2), 
  truncation = trunc_val,
  global_rescale = TRUE
)

morb_tx <- morb_tx_int %>%
  left_join(morb_grappe, by = "morbidity_grappe") %>%
  rename(grappe = morbidity_grappe) %>%
  pivot_longer(cols = ends_with("_int"),
               names_to = "tx_type",
               names_pattern = "(.*)_int",
               values_to = "tx_int") %>%
  st_drop_geometry() %>%
  dplyr::select(grappe, arm, idw_power, tx_type, tx_int) 

# Generate Treatment Intensity Layer --------------------------------------

if (recalc_idw_layer) {
  source(here("scripts", "24-months", "02ax-tx_int-idw-layer.R"))
}

# Permutation Tests -------------------------------------------------------

## Permutation test: Permuting treatment labels of main trial villages
## The treatment intensity layer changes in each permutation, but the morbidity trial
## data remains the same.

perm_tx_int <- idw_rct_perm(
  dist_data = dist_main_morb,
  n_perms = n_permutations,
  power = c(idw_power, 2),
  truncation = trunc_val,
  data_form = "long",
  measure_at_var = "morbidity_grappe",
  treated_at_var = "main_grappe",
  numerator_var = "n_doses",
  n_cores = n_cores
)

perm_tx <- perm_tx_int %>%
  rename(grappe = morbidity_grappe)

# Save Results ------------------------------------------------------------

saveRDS(morb_tx, here("data", "clean", "morb_tx_int.rds"))
saveRDS(perm_tx, here("data", "output", "perm_tx_int.rds"))
