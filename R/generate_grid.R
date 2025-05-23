#' Generate Spatial Grid
#'
#' This function generates a spatial grid over a geographic extent and assigns grid IDs to points.
#' It also computes cell centroids and, for geographic grids (in degrees or minutes), a mapsheet code.
#' Finally, it summarizes specified columns within each grid cell while preserving additional metadata
#' provided via `extra_cols`.
#'
#' @param data A data frame containing point data with x and y coordinates.
#' @param x_col Character. Column name for x-coordinates (default: "x").
#' @param y_col Character. Column name for y-coordinates (default: "y").
#' @param grid_size Numeric. Size of the grid cells. For geographic data (EPSG:4326) this is in degrees.
#' @param sum_col_range Numeric vector. Range of columns to summarize within each grid cell.
#' @param extra_cols Character vector of additional columns to retain in the output (optional).
#' @param crs_epsg Numeric. EPSG code for the coordinate reference system (default: 4326).
#' @param unit Character. One of "deg", "min", "sec", or "m". For geographic data use "deg" (default).
#'
#' @return A list containing:
#'   - `grid`: terra raster object of the generated grid with unique grid IDs.
#'   - `grid_sf`: sf object of the grid as polygons, with added centroid coordinates and (if applicable) mapsheet codes.
#'   - `block_sp`: Data frame summarizing grid cell contents (including extra_cols if provided).
#'
#' @import sf
#' @import terra
#' @import dplyr
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(
#'   x = runif(100, -10, 10),
#'   y = runif(100, -10, 10),
#'   species1 = rpois(100, 5),
#'   species2 = rpois(100, 3),
#'   recordedBy = sample(LETTERS, 100, replace = TRUE)
#' )
#' grid_result <- generate_grid(data, x_col = "x", y_col = "y",
#'                              grid_size = 1, sum_col_range = 3:4,
#'                              extra_cols = c("recordedBy"))
#' print(grid_result$block_sp)
#' plot(grid_result$grid_sf["grid_id"])
generate_grid <- function(data,
                          x_col = "x",
                          y_col = "y",
                          grid_size = 0.5,
                          sum_col_range = NULL,
                          extra_cols = NULL,
                          crs_epsg = 4326,
                          unit = c("deg", "min", "sec", "m")) {

  # Load required packages
  required_packages <- c("sf", "terra", "dplyr", "tidyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.")
    }
  }

  unit <- match.arg(unit)

  # Validate inputs
  if (!all(c(x_col, y_col) %in% names(data))) {
    stop("Specified x or y columns do not exist in the input data frame.")
  }
  if (is.null(sum_col_range)) {
    stop("Please specify `sum_col_range` for columns to sum.")
  }

  # Convert data to sf
  points_sf <- sf::st_as_sf(data, coords = c(x_col, y_col), crs = crs_epsg)

  # Preserve original coords
  data$orig_x <- data[[x_col]]
  data$orig_y <- data[[y_col]]

  # Species columns to summarize
  sum_cols <- names(data)[sum_col_range]

  # ---------------------------------------------------------------------------
  # CASE 1: grid_size == 0 => No grid, assign site_id per unique (x, y)
  # ---------------------------------------------------------------------------
  if (grid_size == 0) {
    message("grid_size = 0: No grid. Assigning grid_id per unique (x, y).")

    # Distinct (orig_x, orig_y)
    coords_df <- data.frame(
      row_id = seq_len(nrow(data)),
      orig_x = data$orig_x,
      orig_y = data$orig_y
    )

    # Each unique coordinate gets a unique grid_id
    unique_coords <- coords_df %>%
      dplyr::distinct(orig_x, orig_y) %>%
      dplyr::mutate(grid_id = as.character(dplyr::row_number()))

    # Join back to the main data
    coords_df <- dplyr::left_join(coords_df, unique_coords, by = c("orig_x", "orig_y"))
    data$grid_id <- coords_df$grid_id

    # Summarize species columns
    block_sp <- data %>%
      dplyr::group_by(grid_id, dplyr::across(dplyr::all_of(extra_cols))) %>%
      dplyr::summarize(
        dplyr::across(all_of(sum_cols), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    # Compute obs_sum, spp_rich
    block_sp <- block_sp %>%
      dplyr::mutate(
        obs_sum  = rowSums(dplyr::select(., all_of(sum_cols)), na.rm = TRUE),
        spp_rich = rowSums(dplyr::select(., all_of(sum_cols)) > 0, na.rm = TRUE)
      )

    # Attach a representative orig_x, orig_y for each grid_id
    block_sp <- block_sp %>%
      dplyr::left_join(
        data %>%
          dplyr::group_by(grid_id) %>%
          dplyr::slice(1) %>%
          dplyr::select(grid_id, orig_x, orig_y),
        by = "grid_id"
      )

    # Reorder columns: (grid_id, orig_x, orig_y, obs_sum, spp_rich, ...)
    front_cols <- c("grid_id", "orig_x", "orig_y", "obs_sum", "spp_rich")
    block_sp <- dplyr::relocate(block_sp, dplyr::all_of(front_cols))

    return(list(
      grid    = NULL,  # No SpatRaster
      grid_sf = NULL,  # No grid polygons
      block_sp = as.data.frame(block_sp)
    ))
  }

  # ---------------------------------------------------------------------------
  # CASE 2: grid_size > 0 => Create grid, assign IDs, compute obs_sum, etc.
  # ---------------------------------------------------------------------------
  message("Generating grid of size ", grid_size, " ", unit, "...")

  # Expand bounding box slightly
  bounds <- sf::st_bbox(points_sf)
  bounds["xmin"] <- floor(bounds["xmin"] / grid_size) * grid_size - 2 * grid_size
  bounds["ymin"] <- floor(bounds["ymin"] / grid_size) * grid_size - 2 * grid_size
  bounds["xmax"] <- ceiling(bounds["xmax"] / grid_size) * grid_size + 2 * grid_size
  bounds["ymax"] <- ceiling(bounds["ymax"] / grid_size) * grid_size + 2 * grid_size
  bb <- bounds

  # Create grid polygons
  grid_polygons <- sf::st_make_grid(sf::st_as_sfc(bb), cellsize = grid_size, what = "polygons")
  grid_sf <- sf::st_sf(geometry = grid_polygons, crs = crs_epsg) %>%
    dplyr::mutate(
      centroid = sf::st_centroid(geometry),
      centroid_lon = sf::st_coordinates(centroid)[, 1],
      centroid_lat = sf::st_coordinates(centroid)[, 2],
      grid_id = as.character(seq_len(nrow(.)))  # character ID
    )

  # Optionally create mapsheet if unit is deg/min
  if (unit %in% c("deg", "min")) {
    lon_int <- floor(grid_sf$centroid_lon)
    lat_int <- floor(grid_sf$centroid_lat)
    lon_dir <- ifelse(grid_sf$centroid_lon >= 0, "E", "W")
    lat_dir <- ifelse(grid_sf$centroid_lat >= 0, "N", "S")
    grid_sf$mapsheet <- sprintf("%s%03d%s%02dBB", lon_dir, abs(lon_int), lat_dir, abs(lat_int))
  } else {
    grid_sf$mapsheet <- NA
  }

  # Create a raster template
  rast_template <- terra::rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = grid_size,
    crs = paste0("EPSG:", crs_epsg)
  )

  # Join each point to a grid cell
  points_join <- sf::st_join(points_sf, grid_sf["grid_id"], join = sf::st_within)
  data$grid_id <- points_join$grid_id

  # Summarize species columns by grid_id
  block_sp <- data %>%
    dplyr::group_by(grid_id, dplyr::across(dplyr::all_of(extra_cols))) %>%
    dplyr::summarize(
      dplyr::across(all_of(sum_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # Compute obs_sum, spp_rich
  block_sp <- block_sp %>%
    dplyr::mutate(
      obs_sum  = rowSums(dplyr::select(., all_of(sum_cols)), na.rm = TRUE),
      spp_rich = rowSums(dplyr::select(., all_of(sum_cols)) > 0, na.rm = TRUE)
    )

  # Join centroid coords & mapsheet from grid_sf
  centroid_info <- grid_sf %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::select(grid_id, centroid_lon, centroid_lat, mapsheet)

  block_sp <- block_sp %>%
    dplyr::left_join(centroid_info, by = "grid_id")

  # Reorder columns: (grid_id, centroid_lon, centroid_lat, mapsheet, obs_sum, spp_rich, ...)
  front_cols <- c("grid_id", "centroid_lon", "centroid_lat", "mapsheet", "obs_sum", "spp_rich")
  block_sp <- dplyr::relocate(block_sp, dplyr::all_of(front_cols))

  # Attach obs_sum & spp_rich to grid_sf => empty cells => NA
  grid_sf <- grid_sf %>%
    dplyr::left_join(
      block_sp %>% dplyr::select(grid_id, obs_sum, spp_rich),
      by = "grid_id"
    )

  # Rasterize each layer separately
  layer_names <- c("grid_id", "obs_sum", "spp_rich")
  r_list <- list()

  grid_vect <- terra::vect(grid_sf)

  for (layer in layer_names) {
    rlay <- terra::rasterize(
      x = grid_vect,
      y = rast_template,
      field = layer,
      background = NA  # NA for empty cells
    )
    names(rlay) <- layer
    r_list[[layer]] <- rlay
  }

  # Stack the single-layer rasters
  rast_final <- do.call(c, unname(r_list))

  return(list(
    grid    = rast_final,  # 3-layer SpatRaster
    grid_sf = grid_sf,     # polygons with grid_id, centroid_lon, centroid_lat, mapsheet, obs_sum, spp_rich
    block_sp = as.data.frame(block_sp)  # summarized data
  ))
}

