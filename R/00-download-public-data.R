#--------------------------------------------
#
# mordor-spillover/R/00-download-public-data.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Download .rds versions of 
# public datasets from osf.io 
# and save them in the local working repository:
# ~/mordor-spillover/data
#
# data are available here: https://osf.io/bmjd3
#
#------------------------------------

library(here)

source(here("R", "00-config.R"))
options(timeout = 600)

if (!dir.exists(path_data_clean)) {
  dir.create(path_data_clean, recursive = TRUE)
}

# OSF ---------------------------------------------------------------------

# Exposure data

## IDW treatment intensity
### https://osf.io/qj64h
morb_tx_int <- osf_retrieve_file("qj64h") %>%
  osf_download(conflict = "overwrite", path = path_data_clean, progress = TRUE)

## Doses in rings
### https://osf.io/zm4dx
morb_rings <- osf_retrieve_file("zm4dx") %>%
  osf_download(conflict = "overwrite", path = path_data_clean, progress = TRUE)

# Outcome data
### https://osf.io/ybj8x
amr_dat <- osf_retrieve_file("ybj8x") %>%
  osf_download(conflict = "overwrite", path = path_data_clean, progress = TRUE)


# Humanitarian Data Exchange ----------------------------------------------

# Niger admin boundaries
### Map downloaded from https://data.humdata.org/dataset/cod-ab-ner
niger_shp_url <- "https://data.humdata.org/dataset/c0e0998c-b45a-4aea-ac06-c1de1d94e596/resource/b2a4cf8d-da46-4f52-bed0-865160470dac/download/ner_admin_boundaries.shp.zip"

niger_zip <- file.path(path_data_clean, "niger_shapefiles.zip")
niger_dir <- file.path(path_data_clean, "niger_shapefiles")

download.file(niger_shp_url, destfile = niger_zip)

if (!dir.exists(niger_dir)) {
  dir.create(niger_dir, recursive = TRUE)
}

unzip(zipfile = niger_zip, exdir = niger_dir)

# Meta high-resolution settlement layer
### Data downloaded from https://data.humdata.org/dataset/highresolutionpopulationdensitymaps-ner#
hrsl_gen_url <- "https://data.humdata.org/dataset/ab6939a8-2546-48db-836e-644150628a2d/resource/909222ae-b255-4ce6-99d5-f79c34136ebc/download/ner_general_2020_csv.zip"
hrsl_u5_url <- "https://data.humdata.org/dataset/ab6939a8-2546-48db-836e-644150628a2d/resource/16904d18-99ac-4282-8656-37d7bca3df71/download/ner_children_under_five_2020_csv.zip"

hrsl_gen_zip <- file.path(path_data_clean, "niger_popdens_general.zip")
hrsl_u5_zip <- file.path(path_data_clean, "niger_popdens_u5.zip")

download.file(hrsl_gen_url, destfile = hrsl_gen_zip)
download.file(hrsl_u5_url, destfile = hrsl_u5_zip)

unzip(zipfile = hrsl_gen_zip, exdir = path_data_clean)
unzip(zipfile = hrsl_u5_zip, exdir = path_data_clean)

# International boundaries shapefile
### Map downloaded from https://data.humdata.org/dataset/global-international-boundaries-osm
int_bound_url <- "https://data.fieldmaps.io/adm0/osm/intl/adm0_polygons.gpkg.zip"

int_bound_zip <- file.path(path_data_clean, "intl_boundaries.zip")
int_bound_dir <- file.path(path_data_clean, "intl_boundaries")

if (!dir.exists(int_bound_dir)) {
  dir.create(int_bound_dir, recursive = TRUE)
}

download.file(int_bound_url, destfile = int_bound_zip)
unzip(zipfile = int_bound_zip, exdir = int_bound_dir)

# Clean up zip files
file.remove(niger_zip, hrsl_gen_zip, hrsl_u5_zip, int_bound_zip)
