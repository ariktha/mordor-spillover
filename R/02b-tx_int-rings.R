#--------------------------------------------
#
# mordor-spillover/R/02b-tx_int-rings.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Calculate the number of doses distributed in rings around 
# morbidity trial villages in main trial villages.
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

# Load distance matrix

dist_main_morb <- readRDS(here("data", "output", "dist_main_morb.RDS")) 

# Doses within concentric rings ------------------------------------------

primary_rings <- tibble(outer_radius = primary_out_radii) %>% 
  mutate(inner_radius = lag(outer_radius, default = 0))

# Compute ring metrics for morbidity trial village locations
ring_morb <- ring_metrics_func(
  dist_long = dist_main_morb, 
  radii_tibble = primary_rings, 
  measure_at_var = "morbidity_grappe",
  only_doses = FALSE)

morb_tx <- ring_morb %>% 
  rename(grappe = morbidity_grappe) %>%
  left_join(morb_grappe, by = "grappe")

morb_tx_sub <- morb_tx %>%
  dplyr::select(grappe, arm, ring_range_text, outer_radius, inner_radius, azithro_doses, placebo_doses)

# Permutations ------------------------------------------------------------

main_perm <- main_grappe %>%
  dplyr::select(grappe, treat_bin) %>%
  rename(treat_bin_obs = treat_bin) %>%
  distinct()

perm_treat_list <- mclapply(1:n_permutations, function(x) {
  sample(main_perm$treat_bin_obs, replace = FALSE)
}, mc.cores = n_cores)

perm_rings_list <- mclapply(perm_treat_list, function(x) {
  main_perm$treat_bin <- x
  main_perm_dist <- dist_main_morb %>%
    rename(grappe = main_grappe) %>%
    dplyr::select(-treat_bin) %>%
    left_join(main_perm, by = "grappe")
  ring_metrics_func(
    dist_long = main_perm_dist,
    radii_tibble = primary_rings,
    measure_at_var = "morbidity_grappe",
    only_doses = TRUE)
}, mc.cores = n_cores)

perm_rings_df <- bind_rows(perm_rings_list, .id = "perm_id") %>%
  rename(grappe = morbidity_grappe)

# Save data --------------------------------------------------------------

saveRDS(morb_tx_sub, here("data", "clean", "morb_rings.rds"))
saveRDS(moirb_tx, here("data", "output", "morb_rings_full.rds"))
saveRDS(perm_rings_df, here("data", "output", "perm_rings.rds"))

