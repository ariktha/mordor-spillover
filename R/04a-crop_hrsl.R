#--------------------------------------------
#
# mordor-spillover/R/04a-crop_hrsl.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Crop high-resolution settlement layer data to the study area, to reduce
# processing time and memory usage.
#
#------------------------------------

library(here)
library(units)

rm(list = ls())

source(here("R", "00-config.R"))


# Load Niger shapefile and filter study area -----------------------------

# Load Niger map from Humanitarian Data Exchange 
# Map downloaded from https://data.humdata.org/dataset/cod-ab-ner on 9 April 2024

niger_shp <- st_read(here("data", "untouched", "niger_shapefiles", "NER_admbnda_adm2_IGNN_20230720.shp"))

# st_crs(niger_shp)

study_admin2 <- c("Falmey", "Boboye", "Loga")
study_shp <- niger_shp %>% dplyr::filter(ADM2_FR %in% study_admin2)

study_union <- st_union(study_shp)


# Create a buffer shape to crop HRSL --------------------------------------

# Transform to a projected CRS suitable for Niger (UTM Zone 31N)
study_union_projected <- st_transform(study_union, proj_crs)

# Buffer in kilometers (convert to meters)
buffer_dist <- set_units(50, "km")
study_buffer <- st_buffer(study_union_projected, dist = buffer_dist)

# Transform back to global CRS
study_buffer <- st_transform(study_buffer, global_crs)
bbox_crop <- st_bbox(study_buffer)

# ggplot() +
#   geom_sf(data = study_union, fill = NA) +
#   geom_sf(data = study_buffer, fill = NA) +
#   coord_sf(xlim = c(bbox_crop["xmin"], bbox_crop["xmax"]), 
#            ylim = c(bbox_crop["ymin"], bbox_crop["ymax"])) +
#   ggspatial::annotation_scale() + theme_minimal()


# Load and crop HRSL data -------------------------------------------------

# Load Meta high-res population density data for Niger from Humanitarian Data Exchange 
# Data downloaded from https://data.humdata.org/dataset/highresolutionpopulationdensitymaps-ner# on 13 November 2024

pop_u5_raw <- data.table::fread(here("data", "untouched", "niger_popdens", "ner_children_under_five_2020.csv"))
pop_gen_raw <- data.table::fread(here("data", "untouched", "niger_popdens", "ner_general_2020.csv"))

pop_u5 <- subset(pop_u5_raw, longitude >= bbox_crop["xmin"] & longitude <= bbox_crop["xmax"] & 
                   latitude >= bbox_crop["ymin"] & latitude <= bbox_crop["ymax"]) %>%
  rename(pop_u5 = ner_children_under_five_2020)

pop_gen <- subset(pop_gen_raw, longitude >= bbox_crop["xmin"] & longitude <= bbox_crop["xmax"] & 
                    latitude >= bbox_crop["ymin"] & latitude <= bbox_crop["ymax"]) %>%
  rename(pop_gen = ner_general_2020)

pop_all <- left_join(pop_gen, pop_u5, by = c("latitude", "longitude"))


# Save cropped HRSL -------------------------------------------------------

saveRDS(pop_all, here("data", "output", "pop_crop.RDS"))

