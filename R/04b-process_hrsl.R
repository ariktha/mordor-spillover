#--------------------------------------------
#
# mordor-spillover/R/04b-process_hrsl.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Process high-resolution settlement layer data to create buffer rings 
# around morbidity grappes and estimate total and under-5 population in these rings
#
#------------------------------------

rm(list = ls())       

library(here)
library(units)
library(tictoc)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))


# Load data ---------------------------------------------------------------

# Load grappe data
morb_gps <- readRDS(here("data", "clean", "morb_gps.RDS"))
morb_grappe_raw <- readRDS(here("data", "clean", "morb_grappe.RDS"))

# Load Meta high-res population density data for Niger from Humanitarian Data Exchange 
# Data downloaded from https://data.humdata.org/dataset/highresolutionpopulationdensitymaps-ner# on 13 November 2024
## Cropped using 04a-crop_hrsl.R script

pop_raw <- readRDS(here("data", "output", "pop_crop.RDS"))

# Load Niger map from Humanitarian Data Exchange 
# Map downloaded from https://data.humdata.org/dataset/cod-ab-ner on 9 April 2024

niger_shp <- st_read(here("data", "untouched", "niger_shapefiles", "NER_admbnda_adm2_IGNN_20230720.shp")) %>%
  mutate(mordor_states = ifelse(ADM2_FR %in% c("Falmey", "Boboye", "Loga"), TRUE, FALSE))

niger_union <- niger_shp %>% 
  dplyr::filter(mordor_states) %>%
  st_union()


# Prep HRSL data ----------------------------------------------------------

pop_dat <- pop_raw %>%
  mutate(grid_id = row_number()) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = global_crs) %>%
  mutate(latitude = st_coordinates(geometry)[, 2],
         longitude = st_coordinates(geometry)[, 1])


# Prep trial data ---------------------------------------------------------

morb_grappe <- morb_grappe_raw %>% 
  left_join(morb_gps, by = "grappe") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = global_crs)

# Create tibble of rings

radii_tib <- tibble(outer_radius = seq(10, 30, by = 10)) %>% 
  mutate(outer_radius = set_units(outer_radius, "km")) %>%
  mutate(inner_radius = outer_radius - set_units(10, "km"))

# Create buffer geometries

morb_buff <- morb_grappe %>%
  dplyr::select(grappe) %>%
  cross_join(radii_tib) %>%
  st_transform(crs = proj_crs) %>%
  mutate(outer_circle = st_buffer(geometry, dist = outer_radius),
         inner_circle = st_buffer(geometry, dist = inner_radius)) %>%
  st_drop_geometry() %>% rowwise() %>%
  mutate(ring_geom = st_difference(outer_circle, inner_circle)) %>% ungroup() %>%
  dplyr::select(-outer_circle, -inner_circle) %>%
  st_as_sf() %>%
  st_transform(crs = global_crs)

# Clip buffer geometries to study region boundaries

morb_buff_study <- morb_buff %>%
  mutate(ring_geom_study = st_intersection(ring_geom, niger_union)) %>% 
  mutate(ring_area_study = set_units(st_area(ring_geom_study), "km^2")) %>% 
  mutate(ring_area_study = drop_units(ring_area_study), outer_radius = drop_units(outer_radius))

# Calculate area of each buffer ring in km^2 and proportion of area in study region

buff_area <- morb_buff %>% 
  mutate(ring_area = set_units(st_area(ring_geom), "km^2")) %>% 
  st_drop_geometry() %>% 
  dplyr::select(outer_radius, ring_area) %>% 
  mutate(across(everything(), drop_units)) %>%
  group_by(outer_radius) %>% summarise(ring_area = mean(ring_area)) %>% ungroup()

# Join buffer area with study area and calculate proportion of buffer area in study region

morb_buff_prop <- morb_buff_study %>%
  left_join(buff_area, by = "outer_radius") %>%
  mutate(prop_in_study_area = ring_area_study/ring_area) %>%
  st_drop_geometry()


# Filter and join HRSL to trial region ------------------------------------

# Bounding boxes for buffers

morb_bbox <- st_bbox(morb_buff)
morb_bbox_study <- st_bbox(morb_buff_study)

# Subset pop based on overlapping bounding box

pop_filtered <- st_crop(pop_dat, morb_bbox)
pop_filtered_study <- st_crop(pop_dat, morb_bbox_study)

# Spatial join to find population within buffers

tic("Spatial Join")
joined <- st_join(pop_filtered, morb_buff, join = st_within)
toc()

tic("Spatial Join: Study region only")
joined_study <- st_join(pop_filtered_study, morb_buff_study, join = st_within)
toc()

# Group by village and summarize population in rings around each morbidity grappe

aggregated <- joined %>% 
  st_drop_geometry() %>%
  group_by(grappe, outer_radius) %>% 
  summarise(total_gen_pop = sum(pop_gen, na.rm = TRUE),
            total_u5_pop = sum(pop_u5, na.rm = TRUE), .groups = "keep") %>% 
  ungroup() %>% drop_na() %>% mutate(outer_radius = drop_units(outer_radius))

aggregated_study <- joined_study %>% 
  st_drop_geometry() %>%
  group_by(grappe, outer_radius) %>% 
  summarise(total_gen_pop = sum(pop_gen, na.rm = TRUE),
            total_u5_pop = sum(pop_u5, na.rm = TRUE), .groups = "keep") %>% 
  ungroup() %>% drop_na() %>% 
  rename(study_area_total_gen_pop = total_gen_pop, study_area_total_u5_pop = total_u5_pop)

aggregated <- aggregated %>% 
  left_join(aggregated_study, by = c("grappe", "outer_radius")) %>%
  left_join(morb_buff_prop %>% dplyr::select(-ring_geom_study), by = join_by(grappe, outer_radius)) %>% 
  dplyr::select(-inner_radius)


# Save the aggregated data ------------------------------------------------

saveRDS(aggregated, here("data", "output", "hrsl_rings.RDS"))
