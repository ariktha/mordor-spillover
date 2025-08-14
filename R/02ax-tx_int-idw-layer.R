# Setup -------------------------------------------------------------------

library(here)
library(raster)

rm(list = ls())       

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))
source(here("R", "02ax-layer-fns.R"))

# Load data ---------------------------------------------------------------

# Load MORDOR morbidity and main trial data

morb_grappe <- readRDS(here("data", "clean", "morb_grappe.RDS"))
main_grappe <- readRDS(here("data", "clean", "main_grappe.RDS"))

# Load gps data

morb_gps <- readRDS(here("data", "clean", "morb_gps.RDS"))
main_gps <- readRDS(here("data", "clean", "main_gps.RDS"))

# Load Niger map from Humanitarian Data Exchange 
# Map downloaded from https://data.humdata.org/dataset/cod-ab-ner on 9 April 2024

niger_shp <- st_read(here("data", "untouched", "niger_shapefiles", "NER_admbnda_adm2_IGNN_20230720.shp"))


# Prep for distance calculation -------------------------------------------

main_grappe <- main_grappe %>% 
  left_join(main_gps, by = "grappe") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = global_crs, remove = FALSE)

# Define study region as Falmey, Boboye, and Loga districts in Niger
niger_sub <- niger_shp %>%
  dplyr::filter(ADM2_FR %in% c("Falmey", "Boboye", "Loga")) %>%
  st_transform(global_crs)

niger_boxgrid <- st_make_grid(niger_sub, n = c(pred_grid_dim, pred_grid_dim), what = "centers") %>%
  st_as_sf(crs = global_crs) %>% rename(geometry = x) %>%
  mutate(grid_id = row_number(),
         grid_long = st_coordinates(geometry)[, 1],
         grid_lat = st_coordinates(geometry)[, 2])

# Crop grid to study region
niger_boxgrid <- niger_boxgrid %>%
  st_filter(niger_sub) 

# Create distance matrix --------------------------------------------------

# Break into Chunks for Distance Matrix Calculation
chunk_size <- 1000
grid_chunks <- split(niger_boxgrid, ceiling(seq_along(niger_boxgrid$grid_id) / chunk_size))

# Parallelize Distance Matrix Computation
dist_main_grid_list <- mclapply(grid_chunks, function(grid_chunk) {
  compute_distance_matrix(
    main_grappe,
    grid_chunk,
    "grappe",
    "grid_id"
  )
}, mc.cores = n_cores)

# Combine Results
dist_main_grid <- bind_rows(dist_main_grid_list) %>%
  left_join(
    main_grappe %>%
      dplyr::select(grappe, treat_bin, n_doses) %>%
      st_drop_geometry(),
    by = "grappe"
  ) %>% mutate(grid_id = as.integer(grid_id)) %>%
  left_join(
    niger_boxgrid %>% dplyr::select(grid_id, grid_long, grid_lat),
    by = "grid_id"
  )

grid_tx_int_rast <- create_idw_raster(
  mask_sf = niger_sub,
  dist_long = dist_main_grid, 
  measure_at_var = "grid_id", 
  numerator_var = "n_doses",
  power = idw_power, 
  truncation = trunc_val
)

grid_tx_int_df <- as.data.frame(grid_tx_int_rast, xy = TRUE)
  
# Save Results ------------------------------------------------------------

saveRDS(grid_tx_int_rast, here("data", "output", "grid_tx_int_rast.RDS"))
saveRDS(grid_tx_int_df, here("data", "output", "grid_tx_int_df.RDS"))
