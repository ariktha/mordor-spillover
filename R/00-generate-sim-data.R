#--------------------------------------------
#
# mordor-spillover/R/00-download-public-data.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Simulate stand-in datasets for
# public release of code
#
# data are available here: https://osf.io/bmjd3
#
#------------------------------------

library(here)

source(here("R", "00-config.R"))

# Generate random coordinates ---------------------------------------------

set.seed(42)

# Number of random locations to generate
n_morb <- 30
n_main <- 594
n_csi <- 200 # There are 183 CSI locations in the 3 districts but closest is not always in the study region

niger_shp <- st_read(file.path(path_data_clean, "niger_shapefiles", "ner_admin2.shp")) %>%
  mutate(mordor_states = ifelse(adm2_name %in% c("Falmey", "Boboye", "Loga"), TRUE, FALSE)) %>%
  dplyr::filter(mordor_states)

niger_union <- st_union(niger_shp)

# Function to generate random points within niger_union

generate_random_points <- function(n_points, polygon) {
  points <- st_sample(polygon, size = n_points, type = "random")
  points_sf <- st_sf(geometry = points)
  coords <- st_coordinates(points_sf)
  coords_df <- as.data.frame(coords)
  coords_df <- coords_df %>%
    rename(longitude = X, latitude = Y)
  return(coords_df)
}

# Generate and save random points

## Morbidity trial GPS coordinates with grappe names

morb_tx_int <- readRDS(file.path(path_data_clean, "morb_tx_int.rds")) %>%
  dplyr::select(grappe, arm) %>%
  distinct()

morb_coords <- generate_random_points(n_morb, niger_union) %>%
  mutate(grappe = morb_tx_int$grappe)

saveRDS(morb_coords, file = file.path(path_data_clean, "morb_gps.rds"))

## Main trial GPS coordinates
main_coords <- generate_random_points(n_main, niger_union) %>%
  mutate(grappe = paste0("main-trial-grappe-", 1:n_main))

saveRDS(main_coords, file = file.path(path_data_clean, "main_gps.rds"))

## CSI GPS coordinates

csi_coords <- generate_random_points(n_csi, niger_union) %>%
  mutate(id = paste0("csi-", 1:n_csi))

saveRDS(csi_coords, file = file.path(path_data_clean, "csi_gps.rds"))


# Sample number of children -----------------------------------------------

set.seed(123)

# Function to sample number of children from a Normal distribution

mean_children <- 210
sd_children <- 167.8572
min_children <- 6

morb_grappe <- morb_tx_int %>%
  mutate(n_children = round(rnorm(n = n(), mean = mean_children, sd = sd_children))) %>%
  mutate(n_children = ifelse(n_children < 0, min_children, round(n_children)),
         treat_bin = ifelse(arm == "azithro", 1, 0))

saveRDS(morb_grappe, file = file.path(path_data_clean, "morb_grappe.rds"))

main_grappe <- data.frame(
  grappe = paste0("main-trial-grappe-", 1:n_main)
) %>%
  mutate(n_children = round(rnorm(n = n(), mean = mean_children, sd = sd_children))) %>%
  mutate(n_children = ifelse(n_children < 0, min_children, round(n_children))) %>%
  mutate(arm = sample(c("azithro", "placebo"), size = n(), replace = TRUE, prob = c(0.5, 0.5))) %>%
  mutate(treat_bin = ifelse(arm == "azithro", 1, 0),
         n_treated = round(0.98*n_children),
         n_doses = n_treated * 4)

saveRDS(main_grappe, file = file.path(path_data_clean, "main_grappe.rds"))
