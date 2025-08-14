#--------------------------------------------
#
# mordor-spillover/R/05a-process_csi_dist.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Estimate distance from morbidity trial villages to the nearest health post
# (CSI) and create a distance matrix.
#
#------------------------------------

rm(list = ls())       

library(here)
library(units)
library(patchwork)
library(ggpubr)

source(here("R", "00-config.R"))

# Load data ---------------------------------------------------------------

# Load MORDOR morbidity and main trial data

morb_grappe <- readRDS(here("data", "clean", "morb_grappe.RDS"))

morb_gps <- readRDS(here("data", "clean", "morb_gps.RDS")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = global_crs)

# Load data on CSI locations

csi_loc_raw <- data.table::fread(here("data", "untouched", "cartesanitaire_niger_entites.csv"))

# Prep CSI locations ------------------------------------------------------

csi_loc <- csi_loc_raw %>% 
  mutate(lat_lon = map(Coordinates, ~ {
    if (!is.na(.x)) {
      # Remove spaces and extract numbers
      clean_coords <- gsub("\\s+", "", .x)
      matches <- regmatches(clean_coords, gregexpr("-?\\d+\\.\\d+", clean_coords))[[1]]
      if (length(matches) >= 2) {
        as.numeric(matches)
      } else {
        c(NA_real_, NA_real_)
      }
    } else {
      c(NA_real_, NA_real_)
    }
  }),
  latitude = map_dbl(lat_lon, 1),
  longitude = map_dbl(lat_lon, 2)) %>%
  dplyr::select(`Dhis2 Id`, latitude, longitude) %>%
  drop_na() %>%
  rename(dhis2_id = `Dhis2 Id`) %>%
  st_as_sf(coords = c("latitude", "longitude"), crs = 4326) %>%
  mutate(latitude = st_coordinates(.)[, 2],
         longitude = st_coordinates(.)[, 1]) %>%
  dplyr::filter(latitude != 0)

# Transform to global CRS

csi_proj <- st_transform(csi_loc, crs = proj_crs)
morb_proj <- st_transform(morb_gps, crs = proj_crs)

# Create distance matrix --------------------------------------------------

morb_csi <- st_distance(morb_proj, csi_proj, by_element = FALSE) %>%
  units::set_units("km") %>%
  as.matrix() %>%
  as.data.frame() %>%
  `colnames<-`(csi_proj[["dhis2_id"]]) %>%
  mutate(grappe = morb_proj[["grappe"]]) %>%
  pivot_longer(-grappe, names_to = "dhis2_id", values_to = "distance") %>%
  mutate(distance = as.numeric(distance))

morb_nearest_csi <- morb_csi %>%
  group_by(grappe) %>%
  slice_min(distance, with_ties = FALSE) %>%
  ungroup()

morb_grappe <- morb_grappe %>% 
  left_join(morb_nearest_csi, by = "grappe")

morb_csi <- morb_grappe %>%
  dplyr::select(grappe, dhis2_id, distance) %>%
  rename(nearest_csi = dhis2_id, dist_csi = distance)

# Save the final data -----------------------------------------------------

saveRDS(morb_csi, here("data", "output", "morb_csi.RDS"))

