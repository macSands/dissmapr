#' Get Environmental Data
#'
#' This function retrieves environmental data for spatial points, either by downloading data from
#' specified sources (e.g., WorldClim, SoilGrids) or using local raster files. Data is extracted for
#' a buffer region around the input points, and missing values are interpolated using surrounding values.
#'
#' @param data A data frame containing spatial points with x and y coordinates.
#' @param buffer_km Numeric. Buffer distance around the input points in kilometers (default: 10).
#' @param var Character. Variable to retrieve (e.g., "bio" for bioclimatic variables, "elev" for elevation).
#' @param res Numeric. Resolution of the environmental data in arc-minutes (default: 2.5).
#' @param path Character. Directory path to store downloaded data or to load local rasters (default: "data/").
#' @param source Character. Data source: "geodata" for downloading or "local" for local rasters (default: "geodata").
#' @param year Numeric. Year for temporal data (e.g., for population or footprint datasets; required for specific datasets).
#' @param model Character. Climate model for CMIP6 projections (e.g., "ACCESS-CM2").
#' @param ssp Numeric. Shared Socioeconomic Pathway for CMIP6 (e.g., 245 for SSP2-4.5).
#' @param time Character. Time period for CMIP6 projections (e.g., "2041-2060").
#' @param depth Character. Depth range for soil data (e.g., "0-5cm").
#' @param stat Character. Statistic to summarize soil data (e.g., "mean").
#'
#' @return A list containing:
#'   - `env_rast`: Raster of the environmental data for the area of interest.
#'   - `sites_sf`: Spatial points as an sf object.
#'   - `env_df`: Data frame combining input data and extracted environmental variables.
#'
#' @import sf
#' @import terra
#' @import dplyr
#' @import geodata
#' @import zoo
#' @export
#'
#' @examples
#' # Example usage:
#' data = data.frame(site_id = 1:5, x = c(10, 12, 14, 16, 18), y = c(20, 22, 24, 26, 28))
#' env_data = get_enviro_data(data, buffer_km = 5, var = "bio", res = 2.5, path = "data/")
#' plot(env_data$env_rast[[1]])
#' points(env_data$sites_sf)
get_enviro_data <- function(data,
                            buffer_km = 10,
                            source = "geodata",
                            var = "bio",
                            res = 2.5, # 10,5,2.5,0.30
                            path = "data/",
                            year = NULL,
                            depth = NULL,
                            stat = "mean",
                            model = NULL,
                            ssp = NULL,
                            time = NULL,
                            temporal = FALSE,
                            time_col = "year",
                            layer_split_char = NULL,
                            sp_cols = NULL,
                            ext_cols = NULL
) {
  required_packages <- c("terra", "sf", "dplyr", "geodata", "zoo")
  lapply(required_packages, function(pkg)
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required but not installed."))

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }

  is_sf_data <- inherits(data, "sf")
  is_polygons <- FALSE
  is_points   <- FALSE

  if (is_sf_data) {
    geom_type <- unique(sf::st_geometry_type(data))
    is_polygons <- any(grepl("POLYGON", geom_type, ignore.case = TRUE))
    is_points   <- any(grepl("POINT",   geom_type, ignore.case = TRUE))
    if (!is_polygons && !is_points) stop("Unsupported geometry type in 'data'.")
    sites_sf <- data
  } else {
    x_cols <- c("x", "lon", "longitude", "decimalLongitude","centroid_lon")
    y_cols <- c("y", "lat", "latitude", "decimalLatitude","centroid_lat")
    find_col <- function(cols) intersect(tolower(names(data)), tolower(cols))[1]
    x <- find_col(x_cols)
    y <- find_col(y_cols)
    if (is.na(x) || is.na(y)) stop("Coordinate columns not found.")
    sites_sf <- sf::st_as_sf(data, coords = c(x, y), crs = 4326, remove = FALSE)
    is_points <- TRUE
  }

  aoi_extent <- if (is_polygons) sf::st_union(sites_sf) else sf::st_buffer(sf::st_convex_hull(sf::st_union(sites_sf)), dist = buffer_km * 1000)

  env_rast <- switch(source,
                     geodata = geodata::worldclim_global(var, res, path),
                     local = terra::rast(var),
                     stop("Unsupported source."))

  env_rast <- terra::crop(env_rast, terra::vect(aoi_extent))
  names(env_rast) <- sub(".*_", "", names(env_rast))

  if (is_polygons) {
    env_data <- terra::extract(env_rast, terra::vect(sites_sf), fun = mean, na.rm = TRUE, ID = TRUE)
    env_df   <- dplyr::bind_cols(as.data.frame(sites_sf), env_data[, -1])
  } else {
    env_data <- terra::extract(env_rast, terra::vect(sites_sf))
    env_df   <- dplyr::bind_cols(as.data.frame(sites_sf), env_data)
  }

  if (!is.null(ext_cols)) env_df <- dplyr::bind_cols(env_df, data[, ext_cols, drop=FALSE])

  env_df <- env_df %>%
    dplyr::mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = FALSE, rule = 2)))

  if (!is.null(sp_cols)) {
    env_df   <- dplyr::select(env_df,  -dplyr::any_of(sp_cols))
    sites_sf <- dplyr::select(sites_sf, -dplyr::any_of(sp_cols))
  }

  list(env_rast = env_rast, sites_sf = sites_sf, env_df = env_df)
}
