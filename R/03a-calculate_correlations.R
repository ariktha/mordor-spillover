#--------------------------------------------
#
# mordor-spillover/R/03a-calculate_correlations.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Calculate correlations between azithromycin treatment intensity and AMR resistance, and
# doses within distance radii and AMR; along with permutation tests.
#
#------------------------------------

rm(list = ls())       

library(here)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))


# Load data ---------------------------------------------------------------

# Load tx intensity data 
morb_tx_int_raw <- readRDS(here("data", "final", "morb_tx_int.RDS"))
perm_tx_int_raw <- readRDS(here("data", "final", "perm_tx_int.RDS"))

# Load ring data
morb_rings_raw <- readRDS(here("data", "final", "morb_rings.RDS"))
perm_rings_raw <- readRDS(here("data", "final", "perm_rings.RDS"))

# Load AMR data
amr_dat <- readRDS(here("data", "final", "amr_dat.RDS"))


# Prepare data ------------------------------------------------------------

# Treatment intensity data

morb_tx_int <- morb_tx_int_raw %>% 
  dplyr::select(grappe, arm, idw_power, tx_type, tx_int) %>%
  left_join(amr_dat, by = "grappe")

perm_tx_int <- perm_tx_int_raw %>% 
  left_join(morb_tx_int %>% dplyr::select(grappe, arm) %>% distinct(), by = "grappe") %>%
  left_join(amr_dat, by = "grappe") 

# Ring data

morb_ring <- morb_rings_raw %>% 
  st_drop_geometry() %>%
  dplyr::select(grappe, arm, ring_range_text, azithro_doses, placebo_doses) %>%
  pivot_longer(cols = c(azithro_doses, placebo_doses), 
               names_to = "tx_type",
               names_pattern = "(.+)_doses",
               values_to = "n_doses") %>%
  left_join(amr_dat, by = "grappe")

perm_ring <- perm_rings_raw %>%
  dplyr::select(grappe, perm_id, ring_range_text, azithro_doses, placebo_doses) %>%
  pivot_longer(cols = c(azithro_doses, placebo_doses), 
               names_to = "tx_type", 
               values_to = "n_doses",
               names_pattern = "(.+)_doses") %>%
  left_join(morb_ring %>% dplyr::select(grappe, arm) %>% distinct(), by = "grappe") %>%
  left_join(amr_dat, by = "grappe", relationship = "many-to-many")


# Treatment intensity and AMR ---------------------------------------------

# Estimate the observed correlations

cor_tx_int <- morb_tx_int %>%
  group_by(arm, phase, tx_type, idw_power, ab_class) %>%
  summarise(cor = safe_correlation(x = tx_int, y = avg_res, est_only = FALSE), .groups = "keep") %>%
  unnest(cor) %>%
  ungroup()

# Set up nested data for parallel processing of permutations

cor_perm_tx_int_nested <- perm_tx_int %>%
  # dplyr::filter(perm_id %in% paste0("perm_", 1:2)) %>%  
  group_by(perm_id, arm, phase, tx_type, idw_power, ab_class) %>%
  nest() %>%
  ungroup()

# Estimate correlations for each permutation

start_time <- Sys.time()

cor_perm_tx_int_nested$cor <- mclapply(
  cor_perm_tx_int_nested$data, 
  function(x) safe_correlation(x$tx_int, x$avg_res, est_only = TRUE), 
  mc.cores = n_cores
)

end_time <- Sys.time()
cat("Time taken for tx intensity correlation calculations:", end_time - start_time, "\n")

# Create a lookup table for observed values 

tx_int_obs_lookup <- cor_tx_int %>%
  dplyr::select(arm, phase, tx_type, idw_power, ab_class, rho_est, R_est) %>%
  {setNames(split(., seq(nrow(.))), 
            paste(.$arm, .$phase, .$tx_type, .$idw_power, .$ab_class, sep = "_"))}

# Process permutations and calculate p-values

tx_int_perm_dists <- cor_perm_tx_int_nested %>%
  dplyr::select(-data) %>%
  unnest(cor) %>%
  mutate(lookup_key = paste(arm, phase, tx_type, idw_power, ab_class, sep = "_")) %>%
  group_split(lookup_key, .keep = TRUE) 

start_time <- Sys.time()

tx_int_pvals <- mclapply(tx_int_perm_dists, function(group_data) {
  tryCatch({
    key <- unique(group_data$lookup_key)
    obs_data <- tx_int_obs_lookup[[key]]
    obs_rho <- obs_data$rho_est[1]
    obs_R <- obs_data$R_est[1]
    
    # Calculate p-values with error handling for each
    rho_pval <- tryCatch({
      get_two_sided_pval(obs = obs_rho, empirical_dist = group_data$rho_est)
    }, error = function(e) NA_real_)
    
    R_pval <- tryCatch({
      get_two_sided_pval(obs = obs_R, empirical_dist = group_data$R_est)
    }, error = function(e) NA_real_)
    
    rho_pval_old <- tryCatch({
      get_two_sided_pval_old(obs = obs_rho, empirical_dist = group_data$rho_est)
    }, error = function(e) NA_real_)
    
    R_pval_old <- tryCatch({
      get_two_sided_pval_old(obs = obs_R, empirical_dist = group_data$R_est)
    }, error = function(e) NA_real_)
    
    tibble(
      arm = obs_data$arm[1],
      phase = obs_data$phase[1], 
      tx_type = obs_data$tx_type[1],
      idw_power = obs_data$idw_power[1],
      ab_class = obs_data$ab_class[1],
      obs_rho = obs_rho,
      obs_R = obs_R,
      perm_rho = list(group_data$rho_est),
      perm_R = list(group_data$R_est),
      rho_perm_p = rho_pval,
      R_perm_p = R_pval
    )
  }, error = function(e) {
    warning(paste("Error processing group:", e$message))
    # Return a row with NAs even when there's an error
    tibble(
      arm = group_data$arm[1],
      phase = group_data$phase[1], 
      tx_type = group_data$tx_type[1],
      idw_power = group_data$idw_power[1],
      ab_class = group_data$ab_class[1],
      obs_rho = obs_rho,
      obs_R = obs_R,
      perm_rho = list(group_data$rho_est),
      perm_R = list(group_data$R_est),
      rho_perm_p = NA_real_,
      R_perm_p = NA_real_
    )
  })
}, mc.cores = n_cores) %>%
  bind_rows()

end_time <- Sys.time()
cat("Time taken for tx intensity p-value calculations:", end_time - start_time, "\n")

# Add observed values to the results

cor_tx_int <- cor_tx_int %>%
  rename(obs_rho = rho_est, obs_R = R_est) %>%
  left_join(tx_int_pvals %>% 
              dplyr::select(arm, phase, tx_type, idw_power, ab_class, obs_rho, obs_R, 
                            rho_perm_p, R_perm_p), 
            by = c("arm", "phase", "tx_type", "idw_power", "ab_class", "obs_rho", "obs_R"))


# Doses in rings and AMR --------------------------------------------------

# Estimate the observed correlations

cor_ring <- morb_ring %>%
  group_by(arm, phase, tx_type, ring_range_text, ab_class) %>%
  summarise(cor = safe_correlation(x = n_doses, y = avg_res, est_only = FALSE), .groups = "keep") %>%
  ungroup() %>%
  unnest(cor)

# Create a lookup table for observed values 

ring_obs_lookup <- cor_ring %>%
  dplyr::select(arm, phase, tx_type, ring_range_text, ab_class, rho_est, R_est) %>%
  {setNames(split(., seq(nrow(.))), 
            paste(.$arm, .$phase, .$tx_type, .$ring_range_text, .$ab_class, sep = "_"))}

# Set up nested data for parallel processing of permutations

cor_perm_ring_nested <- perm_ring %>%
  # dplyr::filter(perm_id %in% c(1:2)) %>%
  group_by(perm_id, arm, phase, tx_type, ring_range_text, ab_class) %>%
  nest() %>%
  ungroup()

start_time <- Sys.time()

cor_perm_ring_nested$cor <- mclapply(
  cor_perm_ring_nested$data, 
  function(x) safe_correlation(x$n_doses, x$avg_res, est_only = TRUE), 
  mc.cores = n_cores
)

end_time <- Sys.time()

cat("Time taken for ring correlation calculations:", end_time - start_time, "\n")

# Process permutations

ring_perm_dists <- cor_perm_ring_nested %>%
  dplyr::select(-data) %>%
  unnest(cor) %>%
  mutate(lookup_key = paste(arm, phase, tx_type, ring_range_text, ab_class, sep = "_")) %>%
  group_split(lookup_key, .keep = TRUE) 

# Calculate p-values for each permutation group

start_time <- Sys.time()

ring_pvals <- mclapply(ring_perm_dists, function(group_data) {
  tryCatch({
    key <- unique(group_data$lookup_key)
    obs_data <- ring_obs_lookup[[key]]
    obs_rho <- obs_data$rho_est[1]
    obs_R <- obs_data$R_est[1]
    
    # Calculate p-values with error handling for each
    rho_pval <- tryCatch({
      get_two_sided_pval(obs = obs_rho, empirical_dist = group_data$rho_est)
    }, error = function(e) NA_real_)
    
    R_pval <- tryCatch({
      get_two_sided_pval(obs = obs_R, empirical_dist = group_data$R_est)
    }, error = function(e) NA_real_)
    
    tibble(
      arm = obs_data$arm[1],
      phase = obs_data$phase[1], 
      tx_type = obs_data$tx_type[1],
      ring_range_text = obs_data$ring_range_text[1],
      ab_class = obs_data$ab_class[1],
      obs_rho = obs_rho,
      obs_R = obs_R,
      perm_rho = list(group_data$rho_est),
      perm_R = list(group_data$R_est),
      rho_perm_p = rho_pval,
      R_perm_p = R_pval
    )
  }, error = function(e) {
    warning(paste("Error processing group:", e$message))
    # Return a row with NAs even when there's an error
    tibble(
      arm = group_data$arm[1],
      phase = group_data$phase[1], 
      tx_type = group_data$tx_type[1],
      ring_range_text = group_data$ring_range_text[1],
      ab_class = group_data$ab_class[1],
      obs_rho = obs_rho,
      obs_R = obs_R,
      perm_rho = list(group_data$rho_est),
      perm_R = list(group_data$R_est),
      rho_perm_p = NA_real_,
      R_perm_p = NA_real_
    )
  })
}, mc.cores = n_cores) %>%
  bind_rows()

end_time <- Sys.time()

cat("Time taken for ring p-value calculations:", end_time - start_time, "\n")

# Add observed values to the results

cor_ring <- cor_ring %>%
  rename(obs_rho = rho_est, obs_R = R_est) %>%
  left_join(ring_pvals %>% 
              dplyr::select(arm, phase, tx_type, ring_range_text, ab_class, obs_rho, obs_R, 
                            rho_perm_p, R_perm_p), 
            by = c("arm", "phase", "tx_type", "ring_range_text", "ab_class", "obs_rho", "obs_R"))


# Save data ---------------------------------------------------------------

## Treatment intensity and AMR correlations

# cor_tx_int has the observed correlations and observational + permutation p-values
# tx_int_perm_dists has the full nested df for each permutation
# tx_int_pvals contains the permutation distributions of correlations along with the the p-values

saveRDS(cor_tx_int, here("data", "final", "cor_tx_int.RDS"))
saveRDS(tx_int_perm_dists, here("data", "detailed", "tx_int_perm_dists.RDS"))
saveRDS(tx_int_pvals, here("data", "final", "tx_int_pvals.RDS"))

## Doses in rings and AMR correlations

saveRDS(cor_ring, here("data", "final", "cor_ring.RDS"))
saveRDS(ring_perm_dists, here("data", "detailed", "ring_perm_dists.RDS"))
saveRDS(ring_pvals, here("data", "final", "ring_pvals.RDS"))
