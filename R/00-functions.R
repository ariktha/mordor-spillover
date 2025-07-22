#--------------------------------------------
#
# mordor-spillover/R/00-functions.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Functions file
#
#------------------------------------


# Helper functions --------------------------------------------------------

# as_km()
#
# Converts a numeric distance from meters to kilometers.
#
# Arguments:
# - x: Numeric. Distance in meters.
#
# Returns:
# - Numeric. Distance in kilometers.

as_km <- function(x) x/1000

# se()
#
# Calculates the standard error of a numeric vector.
#
# Arguments:
# - x: Numeric vector.
#
# Returns:
# - Numeric. Standard error of the mean.

se <- function(x) sqrt(var(x)/length(x))

# safe_correlation()
#
# Computes Spearman and Pearson correlations between two numeric vectors.
#
# Arguments:
# - x, y: Numeric vectors of equal length.
# - est_only: Logical. If TRUE, returns only point estimates; if FALSE, includes p-values.
#
# Details:
# - Uses `cor()` or `cor.test()` depending on `est_only`.
# - Handles errors by returning `NA` values instead of failing.
#
# Returns:
# - A tibble with columns:
#     - rho_est: Spearman correlation.
#     - R_est: Pearson correlation.
#     - (If est_only = FALSE) rho_obs_pval and R_obs_pval: Corresponding p-values.

safe_correlation <- function(x, y, est_only = TRUE) {
  
  if(est_only) {
    
    # Return only the estimate
    tryCatch({
      rho.cor <- cor(x, y, method = "spearman")
      R.cor <- cor(x, y, method = "pearson")
      return(tibble(rho_est = rho.cor, R_est = R.cor))
    }, error = function(e) {
      return(tibble(rho_est = NA_real_, R_est = NA_real_))
    })
    
  } else {
    
    # Return full correlation test results
    tryCatch({
      rho.cor <- cor.test(x, y, method = "spearman", exact = FALSE) # exact = FALSE to avoid warnings about ties
      R.cor <- cor.test(x, y, method = "pearson")
      return(tibble(rho_est = unname(rho.cor$estimate), 
                     rho_obs_pval = rho.cor$p.value,
                     R_est = unname(R.cor$estimate),
                     R_obs_pval = R.cor$p.value))
    }, error = function(e) {
      return(tibble(
        rho_est = NA_real_,
        rho_obs_pval = NA_real_,
        R_est = NA_real_,
        R_obs_pval = NA_real_
      ))
    })
  }
}


# get_two_sided_pval()
#
# Calculates a two-sided p-value from an empirical distribution.
# This function determines how extreme an observed value is relative to
# a null distribution created by permutation testing.
#
# Arguments:
# - obs: Numeric. The observed test statistic.
# - empirical_dist: Numeric vector. The empirical null distribution of the test statistic.
# - plot: Logical. If TRUE, generates a density plot of the empirical distribution
#   with the observed value and median marked. Defaults to FALSE.
#
# Details:
# - For two-sided testing, the function finds the equivalent percentile in the
#   opposite tail of the distribution to ensure proper two-sided p-value calculation.
# - The approach is symmetric around the median of the empirical distribution, not zero.
#
# Requirements:
# - For plotting: `ggplot2` package
#
# Returns:
# - Numeric. The two-sided p-value representing the probability of observing a value
#   as or more extreme than 'obs' under the null distribution.
# - NA if inputs are invalid or insufficient.

get_two_sided_pval <- function(obs, empirical_dist, plot = FALSE) {
  # Check inputs
  if (!is.numeric(obs) || length(obs) != 1) {
    warning("Observed value (obs) must be a single numeric value")
    return(NA_real_)  # Return 1 (conservative) for invalid inputs
  }
  
  if(is.na(obs) || length(empirical_dist) == 0) {
    return(NA_real_)
  }
  
  if (!is.numeric(empirical_dist) || length(empirical_dist) < 2) {
    warning("Empirical distribution must be a numeric vector with sufficient samples")
    return(NA_real_)
  }
  
  # Remove NA, NaN, and Inf values
  empirical_dist <- empirical_dist[is.finite(empirical_dist)]
  
  # Check if we still have enough values
  if (length(empirical_dist) < 2) {
    warning("Not enough finite values in empirical distribution")
    return(NA_real_)
  }
  
  # Check if there's any variation in the empirical distribution
  if (var(empirical_dist) < .Machine$double.eps) {
    warning("No variation in empirical distribution")
    # If obs matches the constant value, return 1, otherwise 0
    return(ifelse(abs(obs - empirical_dist[1]) < .Machine$double.eps, 1, 0))
  }
  
  N <- length(empirical_dist)
  median_val <- median(empirical_dist)
  
  obs_distance <- abs(obs - median_val)
  empirical_distances <- abs(empirical_dist - median_val)
  
  # Count how many empirical values are as or more extreme
  p_value <- sum(empirical_distances >= obs_distance) / N
  
  # Ensure the p-value does not exceed 1
  p_value <- min(p_value, 1)
  
  # Check for NaN
  if (is.na(p_value)) {
    warning("P-value calculation resulted in NA")
    return(NA_real_)  # Return NA
  }
  
  # Create plot if requested
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 package is required for plotting but not available. Skipping plot.")
    } else {
      plot <- ggplot2::ggplot() + 
        ggplot2::geom_density(ggplot2::aes(x = empirical_dist)) + 
        ggplot2::geom_vline(xintercept = median_val) + 
        ggplot2::geom_vline(xintercept = obs, color = "blue", linetype = "dashed") + 
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Empirical Distribution with Observed Value",
          x = "Value",
          y = "Density",
          caption = paste("Two-sided p-value:", round(p_value, 4))
        )
      
      print(plot)
    }
  }
  
  return(p_value)
}


# compute_distance_matrix()
#
# Computes a pairwise distance matrix between two spatial datasets.
# Projects both datasets to a specified CRS, calculates distances (in kilometers), 
# reshapes to long format, and returns a tidy data frame.
#
# Arguments:
# - df1: First spatial data frame (sf object) with geometry.
# - df2: Second spatial data frame (sf object) with geometry.
# - colname1: Column name (string) for the identifier in df1.
# - colname2: Column name (string) for the identifier in df2.
#
# Requirements:
# - df1 and df2 must be valid `sf` spatial objects.
# - External objects: proj_crs: A CRS to which both datasets are projected for distance measurements in meters.
# - Packages: `sf`, `dplyr`, `tidyr`, `units`
#
# Returns:
# - A tidy data frame with three columns:
#     - Identifier from df1 (colname1)
#     - Identifier from df2 (colname2)
#     - Distance in kilometers (numeric)

compute_distance_matrix <- function(df1, df2, colname1, colname2) {
  df1_proj <- st_transform(df1, crs = proj_crs)
  df2_proj <- st_transform(df2, crs = proj_crs)
  
  # Precompute identifiers
  id1 <- df1[[colname1]]
  id2 <- df2[[colname2]]
  
  distance_matrix <- st_distance(df1_proj, df2_proj, by_element = FALSE) %>%
    units::set_units("km") %>%
    as.data.frame()
  
  colnames(distance_matrix) <- id2
  distance_matrix[[colname1]] <- id1
  
  distance_matrix %>%
    pivot_longer(
      cols = -!!sym(colname1),
      names_to = colname2,
      values_to = "distance"
    ) %>%
    mutate(distance = as.numeric(distance))
}

# calculate_idw()
#
# Calculates inverse distance weighted (IDW) treatment intensities from a pairwise distance matrix.
# Supports multiple power values and returns results in long format with one row per location per power.
#
# Arguments:
# - dist_long: A data frame of pairwise distances in long format, containing at least:
#     - 'distance': Numeric distance values.
#     - 'treat_bin': Binary treatment indicator (1 = azithromycin, 0 = placebo).
#     - [numerator_var]: Column to use as the numerator in weighting (e.g., 'n_doses').
# - measure_at_var: String. Column name used to group results (e.g., location ID).
# - numerator_var: String. Column name for the numerator used in IDW weighting.
# - power: Numeric scalar or vector. One or more power values used in the IDW formula.
# - truncation: Numeric vector of length 2:
#     - truncation[1]: Minimum distance threshold.
#     - truncation[2]: Replacement value for distances below this threshold.
#     (Defaults to no truncation: c(0, 0))
# - global_rescale: Logical. If TRUE, rescales intensities globally; if FALSE, rescales within each group.
#
# Details:
# - Computes weights as 1 / distance^power (or using a truncated distance if specified).
# - Aggregates weighted intensities for azithromycin and placebo treatment arms.
# - Rescales intensity values using `scales::rescale()` for comparability across groups or globally.
# - If multiple power values are supplied, computes and combines results across all values.
#
# Requirements:
# - Packages: `dplyr`, `scales`, `purrr`
#
# Returns:
# - A data frame in long format, grouped by `measure_at_var` and `power`, containing:
#     - azithro_int: Rescaled IDW intensity for azithromycin.
#     - placebo_int: Rescaled IDW intensity for placebo.
#     - power: The power value used in the IDW computation.
#
# Notes:
# - Designed to integrate with permutation functions such as `idw_rct_perm()`.
# - Automatically filters out NA and infinite values before summarizing results.


calculate_idw <- function(dist_long, measure_at_var, numerator_var, power, 
                          truncation = c(0, 0), global_rescale = TRUE) {
  # Input validation
  required_cols <- c("distance", numerator_var, "treat_bin", measure_at_var)
  missing_cols <- setdiff(required_cols, names(dist_long))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (any(power <= 0)) {
    stop("Power parameter(s) must be positive")
  }
  
  if (length(truncation) != 2) {
    stop("Truncation must be a numeric vector of length 2")
  }
  
  # Wrapper function to apply calculation at a single power
  calc_single_power <- function(p) {
    if (truncation[1] == 0 && truncation[2] == 0) {
      weighted <- dist_long %>%
        dplyr::mutate(
          weight = 1 / (.data$distance^p),
          weighted_azithro = .data[[numerator_var]] * .data$treat_bin * weight,
          weighted_placebo = .data[[numerator_var]] * (1 - .data$treat_bin) * weight
        )
    } else {
      weighted <- dist_long %>%
        dplyr::mutate(
          dist_idw = dplyr::if_else(.data$distance < truncation[1], truncation[2], .data$distance),
          weight = 1 / (dist_idw^p),
          weighted_azithro = .data[[numerator_var]] * .data$treat_bin * weight,
          weighted_placebo = .data[[numerator_var]] * (1 - .data$treat_bin) * weight
        )
    }
    
    # Filter out problematic values
    weighted <- weighted %>%
      dplyr::filter(
        is.finite(weighted_azithro),
        is.finite(weighted_placebo),
        !is.na(weighted_azithro),
        !is.na(weighted_placebo)
      )
    
    # Summarize and rescale
    summarized <- weighted %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(measure_at_var))) %>%
      dplyr::summarize(
        azithro_int = sum(weighted_azithro),
        placebo_int = sum(weighted_placebo),
        .groups = "drop"
      )
    
    if (global_rescale) {
      summarized <- summarized %>%
        dplyr::mutate(dplyr::across(dplyr::ends_with("_int"), scales::rescale))
    } else {
      summarized <- summarized %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(measure_at_var))) %>%
        dplyr::mutate(dplyr::across(dplyr::ends_with("_int"), scales::rescale)) %>%
        dplyr::ungroup()
    }
    
    summarized %>% dplyr::mutate(idw_power = p)
  }
  
  # Handle single or multiple power values
  result <- purrr::map_dfr(power, calc_single_power)
  
  return(result)
}

# Permutation tests -------------------------------------------------------

# idw_rct_perm()
#
# Performs permutation testing for IDW-based treatment intensity metrics,
# allowing for one or more values of the power parameter.
# For each permutation, treatment labels are shuffled and IDW intensities
# are recalculated. If multiple power values are provided, the function
# computes results for each power and returns results in long format.
#
# Arguments:
# - dist_data: Long-format distance data frame with treatment assignments.
# - n_perms: Integer. Number of permutations to generate.
# - power: Numeric scalar or vector. Power parameter(s) for inverse distance weighting.
# - truncation: Numeric vector of length 2 for distance truncation (min threshold, replacement).
# - data_form: String. Output format; "long" (default) reshapes by treatment type.
# - measure_at_var: String. Variable for location at which intensity is measured.
# - treated_at_var: String. Variable for location receiving treatment.
# - numerator_var: String. Variable to use for weighting numerator (e.g., "n_doses").
# - n_cores: Integer. Number of cores for parallel computation.
#
# Details:
# - For each permutation and power, treatment labels are randomly reassigned
#   and IDW values recalculated using `calculate_idw()`.
# - If `power` is a vector, results are computed for all values and returned stacked.
#
# Requirements:
# - Depends on: `calculate_idw()`
# - Packages: `dplyr`, `purrr`, `parallel`, `tidyr`
#
# Returns:
# - A tidy data frame of permutation results, with one row per permutation,
#   per power value, and optionally per treatment type (if `data_form = "long"`).
#   Includes columns:
#     - perm_id: Permutation index
#     - power: Power value used for IDW
#     - tx_type: Treatment type ("azithro", "placebo") if in long format
#     - tx_int: Rescaled treatment intensity
#     - [measure_at_var]: Grouping location


idw_rct_perm <- function(dist_data, n_perms, power, truncation, data_form = "long",
                         measure_at_var = "morbidity_grappe", treated_at_var = "main_grappe",
                         numerator_var = "n_doses", n_cores = 1) {
  
  treat_sym <- sym("treat_bin")
  treated_sym <- sym(treated_at_var)
  
  set.seed(123)  # For reproducibility
  
  # Generate permuted treatment labels once, reused across power values
  perm_treat_list <- mclapply(1:n_perms, function(x) {
    sample(dist_data$treat_bin, replace = FALSE)
  }, mc.cores = n_cores)
  
  # Loop over power values
  all_results <- purrr::map_dfr(power, function(p) {
    # Apply permutations at each power level
    perm_results <- mclapply(perm_treat_list, function(treat_labels) {
      dist_data$treat_bin <- treat_labels
      calculate_idw(
        dist_data,
        measure_at_var = measure_at_var,
        numerator_var = numerator_var,
        power = p,
        truncation = truncation
      )
    }, mc.cores = n_cores)
    
    names(perm_results) <- paste0("perm_", seq_along(perm_results))
    
    map_dfr(perm_results, as_tibble, .id = "perm_id") %>%
      mutate(idw_power = p)
  })
  
  # Reshape if requested
  if (data_form == "long") {
    all_results <- all_results %>%
      pivot_longer(
        cols = ends_with("_int"),
        names_to = "tx_type",
        values_to = "tx_int"
      ) %>%
      mutate(tx_type = ifelse(tx_type == "azithro_int", "azithro", "placebo"))
  }
  
  return(all_results)
}


# permutation_rho <- function(perm_dat, obs_dat, perm_type = "rct") {
#   
#   obs_corr <- obs_dat %>% 
#     group_by(arm, phase, tx_int_type) %>% 
#     summarise(cor = unname(cor.test(avg_res, tx_int, method = "spearman")$estimate), 
#               .groups = "drop")
#   
#   corr_perm_dat <- perm_dat %>% 
#     group_by(perm_id, arm, phase, tx_int_type) %>% 
#     summarize(perm_cor = unname(cor.test(avg_res, tx_int, method = "spearman")$estimate), 
#               .groups = "drop")
#   
#   corr_perm <- corr_perm_dat %>%
#     left_join(obs_corr, by = c("arm", "phase", "tx_int_type")) %>% 
#     group_by(arm, phase, tx_int_type) %>% 
#     summarise(perm_pval = get_two_sided_pval_new(unique(cor), perm_cor), 
#               .groups = "drop") %>% 
#     left_join(obs_corr, by = c("arm", "phase", "tx_int_type")) %>%
#     mutate(perm_type = perm_type)
#   
#   return(corr_perm)
# }


# Compute ring metrics ----------------------------------------------------

# ring_metrics_func()
#
# Computes treatment metrics (e.g., doses, children treated) within concentric distance bands (rings).
#
# Arguments:
# - dist_long: Long-format distance matrix with treatment and outcome variables.
# - radii_tibble: Tibble with columns `inner_radius` and `outer_radius`.
# - measure_at_var: String. Variable name for where metrics are measured.
# - only_doses: Logical. If TRUE, returns only doses by treatment arm.
#
# Requirements:
# - Depends on: `as_km()` indirectly (via distance interpretation)
# - Packages: `dplyr`, `sf`, `tidyr`
#
# Returns:
# - A data frame with aggregated metrics for each ring per measurement unit.
# - Includes ring metadata (radius bounds and text labels).


ring_metrics_func <- function(dist_long, radii_tibble, measure_at_var = "grappe", only_doses = FALSE) {
  
  # Validate radii_tibble input
  if (!all(c("inner_radius", "outer_radius") %in% colnames(radii_tibble))) {
    stop("radii_tibble must contain 'inner_radius' and 'outer_radius' columns.")
  }
  
  radii_tibble <- radii_tibble %>%
    dplyr::filter(inner_radius >= 0 & outer_radius > 0) %>% # Inner and outer radii must be positive
    dplyr::filter(inner_radius < outer_radius) %>%         # Inner radius must be less than outer radius
    mutate(ring_rank = row_number())                       # Assign a unique rank to each ring
  
    dist_long <- st_drop_geometry(dist_long) %>%
    cross_join(radii_tibble)

  
  if(only_doses) {
    
    ring_tbl <- dist_long %>%
      mutate(
        placebo_bin = ifelse(treat_bin == 0, 1, 0), 
        dist_bin = ifelse(distance > inner_radius & distance <= outer_radius, 1, 0),
        placebo_doses = ifelse(dist_bin == 1, n_doses * placebo_bin, 0),
        azithro_doses = ifelse(dist_bin == 1, n_doses * treat_bin, 0),
      ) %>%
      group_by(!!sym(measure_at_var), ring_rank) %>%
      summarise(
        placebo_doses = sum(placebo_doses),
        azithro_doses = sum(azithro_doses), 
        .groups = "keep"
      ) %>%
      ungroup()
    
  } else {
  
    ring_tbl <- dist_long %>%
      mutate(
        placebo_bin = ifelse(treat_bin == 0, 1, 0), 
        dist_bin = ifelse(distance > inner_radius & distance <= outer_radius, 1, 0),
        placebo_doses = ifelse(dist_bin == 1, n_doses * placebo_bin, 0),
        azithro_doses = ifelse(dist_bin == 1, n_doses * treat_bin, 0),
        # placebo_treated = ifelse(dist_bin == 1, n_treated * placebo_bin, 0),
        # azithro_treated = ifelse(dist_bin == 1, n_treated * treat_bin, 0),
        placebo_children = ifelse(dist_bin == 1, n_children * placebo_bin, 0),
        azithro_children = ifelse(dist_bin == 1, n_children * treat_bin, 0),
        placebo_grappes = ifelse(dist_bin == 1, placebo_bin, 0),
        azithro_grappes = ifelse(dist_bin == 1, treat_bin, 0),
        total_doses = ifelse(dist_bin == 1, n_doses, 0),
        total_grappes = ifelse(dist_bin == 1, 1, 0)
      ) %>%
      group_by(!!sym(measure_at_var), ring_rank) %>%
      summarise(
        placebo_doses = sum(placebo_doses),
        azithro_doses = sum(azithro_doses), 
        placebo_grappes = sum(placebo_grappes),
        azithro_grappes = sum(azithro_grappes), 
        # placebo_treated = sum(placebo_treated),
        # azithro_treated = sum(azithro_treated),
        placebo_children = sum(placebo_children),
        azithro_children = sum(azithro_children),
        total_doses = sum(total_doses),
        total_grappes = sum(total_grappes),
        # total_treated = placebo_treated + azithro_treated,
        total_children = placebo_children + azithro_children,
        .groups = "keep"
      ) %>%
      ungroup()
  }
  
  ring_tbl <- ring_tbl  %>%
    left_join(radii_tibble, by = "ring_rank") %>%
    mutate(
      ring_range_text = paste0(inner_radius, " - ", outer_radius, " km")
    )
}

