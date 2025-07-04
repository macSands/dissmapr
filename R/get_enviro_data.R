#' Retrieve, crop, and link environmental rasters to sampling sites
#'
#' `get_enviro_data()` automates five steps:
#' 1. **AOI build** – buffers the convex hull of all input points.
#' 2. **Raster acquisition** – downloads or loads a multi-layer stack
#'    (WorldClim, SoilGrids, footprint, population, or user-supplied).
#' 3. **Cropping** – trims the stack to the buffered AOI.
#' 4. **Extraction & gap fill** – pulls raster values at each site and
#'    linearly interpolates isolated NAs.
#' 5. **Assembly** – returns a tidy *site × environment* table plus the cropped
#'    raster and an `sf` layer of sites.
#'
#' @section Supported raster sets (`var` when `source = "geodata"`):
#' * `"bio"`  19 WorldClim v2.1 bioclim variables
#' * `"elev"` WorldClim elevation layer
#' * `"footprint"` *Global Human Footprint* (specify `year`)
#' * `"population"` WorldPop population counts (specify `year` & `res`)
#' * `"soil_world"` ISRIC SoilGrids layers (specify `depth`, `stat`)
#'
#' @param data      Data frame of spatial points; must include lon/lat columns
#'                  such as `"x","y"` or `"centroid_lon","centroid_lat"`.
#' @param buffer_km Buffer width (km) for the AOI. Default 10.
#' @param source    `"geodata"` (default) to fetch layers via **geodata** or
#'                  `"local"` to supply a local `SpatRaster` or file path.
#' @param var       Raster product to download (see list above) or ignored when
#'                  `source = "local"`.
#' @param res       Resolution (arc-min) for WorldClim/WorldPop layers.
#' @param path      Cache folder for downloaded rasters (created if absent).
#' @param year, model, ssp, time  Optional arguments for time-stamped products
#'                  (human footprint, population, CMIP6, etc.).
#' @param depth, stat  Arguments passed to `geodata::soil_world()`.
#' @param sp_cols   **Columns to drop** from the final table (e.g. a large
#'                  species matrix). Accepts names or numeric indices *relative
#'                  to `data`*.
#' @param ext_cols  **Columns to append** verbatim (e.g. `"obs_sum","spp_rich"`).
#'
#' @return A list with
#' * `env_rast`  `SpatRaster` – cropped environmental stack (all layers)
#' * `sites_sf`  `sf` POINT layer (WGS-84) of the input sites
#' * `env_df`    Tibble with site ID, coordinates, every raster variable,
#'               plus any columns requested in `ext_cols`
#'
#' @details
#' WorldClim BIO layers arrive as a multi-file *SpatRasterDataset*; the function
#' concatenates all sub-rasters so that **bio01 – bio19** are returned together.
#' Missing cells inside the AOI are filled with linear interpolation
#' (`zoo::na.approx`).
#'
#' @examples
#' # Example usage:
#' data = data.frame(site_id = 1:5, x = c(10, 12, 14, 16, 18), y = c(20, 22, 24, 26, 28))
#' env_data = get_enviro_data(data, buffer_km = 5, var = "bio", res = 2.5, path = "data/")
#' plot(env_data$env_rast[[1]])
#' points(env_data$sites_sf)
#'
#' @import sf terra dplyr geodata zoo
#' @export
get_enviro_data <- function(data,
                            buffer_km = 10,
                            source    = "geodata",
                            var       = "bio",      # "bio", "elev", ...
                            res       = 2.5,
                            path      = "data/",
                            year      = NULL, depth = NULL, stat = "mean",
                            model     = NULL, ssp = NULL, time = NULL,
                            sp_cols   = NULL,       # drop
                            ext_cols  = NULL) {     # append

  ## ── deps ────────────────────────────────────────────────────────────────
  for (pkg in c("terra","sf","dplyr","geodata","zoo"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required but not installed.")

  if (!dir.exists(path)) dir.create(path, recursive = TRUE)

  ## ── identify coord columns ──────────────────────────────────────────────
  x_col <- intersect(tolower(names(data)),
                     c("x","lon","longitude","decimallongitude","centroid_lon"))[1]
  y_col <- intersect(tolower(names(data)),
                     c("y","lat","latitude","decimallatitude","centroid_lat"))[1]
  id_col <- intersect(names(data), c("site_id","grid_id"))[1]

  if (is.na(x_col) || is.na(y_col))
    stop("Coordinate columns not found in `data`.")

  message("• Using coord cols: ", x_col, ", ", y_col)

  ## ── build AOI ───────────────────────────────────────────────────────────
  cols_xy <- c(id_col, x_col, y_col); cols_xy <- cols_xy[!is.na(cols_xy)]
  data_xy <- data |>
    dplyr::select(dplyr::all_of(cols_xy)) |>
    dplyr::distinct()

  sites_sf <- sf::st_as_sf(data_xy, coords = c(x_col, y_col), crs = 4326)
  aoi      <- sf::st_buffer(sf::st_convex_hull(sf::st_union(sites_sf)),
                            buffer_km * 1e3)

  message("• AOI built and buffered by ", buffer_km, " km")

  ## ── download / load rasters ─────────────────────────────────────────────
  message("• Acquiring raster stack")
  env_rast <- switch(source,
                     geodata = switch(var,
                                      bio  = geodata::worldclim_global("bio",  res, path),
                                      elev = geodata::worldclim_global("elev", res, path),
                                      footprint = geodata::footprint(year, path),
                                      population = geodata::population(year, res, path),
                                      soil_world = geodata::soil_world(var, depth, stat, path),
                                      stop("Unsupported `var` for geodata source.")),
                     local = {
                       if (inherits(var, "SpatRaster")) var
                       else if (file.exists(var)) terra::rast(var)
                       else stop("`var` must be a SpatRaster or file path.")},
                     stop("`source` must be 'geodata' or 'local'.")
  )

  ## WorldClim BIO arrives as SpatRasterDataset (3 files * 6–7 layers each)
  if (inherits(env_rast, "SpatRasterDataset")) {
    message("  - Merging ", length(env_rast), " raster files")
    env_rast <- do.call(terra::c, lapply(env_rast, terra::rast))
  }

  message("  - Total layers: ", terra::nlyr(env_rast))

  ## crop and rename
  env_rast <- terra::crop(env_rast, terra::vect(aoi))
  if (var == "bio")
    names(env_rast) <- sprintf("bio%02d", seq_len(terra::nlyr(env_rast)))

  ## ── extract values ──────────────────────────────────────────────────────
  message("• Extracting raster values at ", nrow(sites_sf), " points")
  vals_df <- terra::extract(env_rast, terra::vect(sites_sf), df = TRUE) |>
    as.data.frame()

  message("  - Layers extracted: ", ncol(vals_df) - 1)  # first col is ID

  vals_df <- vals_df[ , -1, drop = FALSE]   # drop ID column
  env_df  <- dplyr::bind_cols(data_xy, vals_df)

  ## ── interpolate small gaps ──────────────────────────────────────────────
  env_df <- dplyr::mutate(
    env_df,
    dplyr::across(where(is.numeric),
                  ~ zoo::na.approx(.x, na.rm = FALSE, rule = 2)))

  ## ── tidy: drop species cols, add extra cols ──────────────────────────────
  if (!is.null(sp_cols)) {
    sp_cols_names <- if (is.numeric(sp_cols)) names(data)[sp_cols] else sp_cols
    env_df <- dplyr::select(env_df, -dplyr::any_of(sp_cols_names))
  }

  if (!is.null(ext_cols))
    env_df <- dplyr::bind_cols(env_df, data[, ext_cols, drop = FALSE])

  message("• Final env_df cols: ", paste(names(env_df), collapse = ", "))

  list(env_rast = env_rast,
       sites_sf = sites_sf,
       env_df   = env_df)
}
