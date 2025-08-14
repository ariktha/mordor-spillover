# create_idw_raster()
#
# Generates a raster surface of inverse distance weighted (IDW) azithromycin treatment intensity
# using a set of spatial points and a mask geometry.
#
# Arguments:
# - mask_sf: An `sf` polygon object used to define the spatial extent and shape of the output raster.
# - dist_long: A long-format data frame containing pairwise distances and treatment variables.
# - measure_at_var: String. Column name identifying the unit at which intensities are measured (e.g., "grid_id").
# - numerator_var: String. Column used as the numerator in the IDW weighting (e.g., "n_doses").
# - power: Numeric scalar or vector. Exponent(s) for inverse distance weighting.
# - truncation: Numeric vector of length 2:
#     - truncation[1]: Minimum distance threshold.
#     - truncation[2]: Replacement value for values below the threshold.
#     (Defaults to c(0, 0) for no truncation.)
# - global_rescale: Logical. If TRUE, rescales intensities globally; if FALSE, rescales within group. Default is TRUE.
#
# Details:
# - First calculates IDW azithromycin intensities using `calculate_idw()`.
# - Joins spatial coordinates (`grid_long`, `grid_lat`) to the IDW values using `grid_id`.
# - Converts the intensity values to a raster layer clipped to the `mask_sf` boundary.
#
# Requirements:
# - Depends on: `calculate_idw()`, `points_to_raster()`
# - Packages: `dplyr`, `sf`
#
# Returns:
# - A raster object representing the spatial distribution of azithromycin IDW intensity.
#
# Notes:
# - Only the azithromycin intensity (`azithro_int`) is used in the raster surface.
# - Assumes `dist_long` includes columns: `grid_id`, `grid_long`, and `grid_lat`.

create_idw_raster <- function(mask_sf, dist_long, measure_at_var, numerator_var, power, 
                              truncation = c(0, 0), global_rescale = TRUE) {
  
  idw_data <- calculate_idw(dist_long = dist_long, 
                            measure_at_var, 
                            numerator_var, 
                            power, 
                            truncation = truncation, 
                            global_rescale = global_rescale)
  
  grid_coords <- dist_long %>% distinct(grid_id, grid_long, grid_lat)
  
  idw_data <- left_join(idw_data, grid_coords, by = c("grid_id" = measure_at_var))
  
  points_to_raster(
    x = idw_data$grid_long,
    y = idw_data$grid_lat,
    z = idw_data$azithro_int,
    mask1 = mask_sf
  )
  
}

#------------------------------------
# Geospatial helper functions from Ben Arnold's repo - https://osf.io/cxb5e 
# Geographic pair matching in large-scale cluster randomized trials
# Functions from file https://osf.io/eqcvh

# points_to_raster()
#
# convert predicted surface of points 
# to a raster to make it a regular
# set of tiles (rather than points), 
# this function takes as input
# an object of class sf dataframe
# in particular, the results of the function
# above, get_grid_preds()
#
# @x,y,z  : longitude, latitude, and some value to plot in raster
# @mask1  : mask the raster by a variable? if NULL, skips.
# @mask2  : mask the raster by a second variable? if NULL skips.
# @crop1  : crop the raster by the first variable?  if NULL, skips

points_to_raster <- function(x, y, z, mask1 = NULL, mask2 = NULL, crop1 = NULL) {
  
  # make raster from points
  pred_raster <- raster::rasterFromXYZ(xyz = data.frame(x = x, y = y, z = z), crs = 4326)
  
  # mask and crop raster if specified
  if(!is.null(mask1)) {
    pred_raster <- mask(x = pred_raster, mask = mask1)
  }
  if(!is.null(mask2)) {
    pred_raster <- mask(x = pred_raster, mask = mask2)
  }
  if(!is.null(crop1)) {
    pred_raster <- crop(x = pred_raster, y = extent(crop1))
  }
  
  return(pred_raster)
}


# raster_to_tibble()
#
#' Transform raster as data.frame to be later used with ggplot
#' Modified from rasterVis::gplot
#' From: https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r
#'
#' @param x A Raster* object
#' @param maxpixels Maximum number of pixels to use
#'
#' @details rasterVis::gplot is nice to plot a raster in a ggplot but
#' if you want to plot different rasters on the same plot, you are stuck.
#' If you want to add other information or transform your raster as a
#' category raster, you can not do it. With `SDMSelect::gplot_data`, you retrieve your
#' raster as a tibble that can be modified as wanted using `dplyr` and
#' then plot in `ggplot` using `geom_tile`.
#' If Raster has levels, they will be joined to the final tibble.
#'
#' @export

raster_to_tibble <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- tibble::as_tibble(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  return(dat)
}


# grid_pred_to_tibble()
#
# Wrapper function to get predictions using get_grid_preds(),
# points_to_raster(), and raster_to_tibble()
#
# the model fit must have variables 
# "lon" and "lat" as the x,y coordinates

grid_pred_to_tibble <- function(model_fit, pred_area_boxgrid){
  grid_pred <- get_grid_preds(input_grid = pred_area_boxgrid, spamm_model_fit = model_fit)
  
  grid_pred_coords <- st_coordinates(grid_pred)
  grid_pred_raster <- points_to_raster(
    x = grid_pred_coords[,1], y = grid_pred_coords[,2], z = grid_pred$pred,
    mask1 = niger_sub, crop1 = niger_sub)
  
  grid_pred_rastib <- raster_to_tibble(grid_pred_raster)
  
  return(grid_pred_rastib)
}

# grid_pred_to_raster()
#
# Wrapper function to get predictions using get_grid_preds(),
# points_to_raster()
#
# the model fit must have variables 
# "lon" and "lat" as the x,y coordinates

grid_pred_to_raster <- function(model_fit, pred_area_boxgrid){
  grid_pred <- get_grid_preds(input_grid = pred_area_boxgrid, spamm_model_fit = model_fit)
  
  grid_pred_coords <- st_coordinates(grid_pred)
  grid_pred_raster <- points_to_raster(
    x = grid_pred_coords[,1], y = grid_pred_coords[,2], z = grid_pred$pred,
    mask1 = niger_sub, crop1 = niger_sub)
  
  return(grid_pred_raster)
}

