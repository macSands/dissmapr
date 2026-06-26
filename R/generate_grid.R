#' Generate Spatial Grid and Gridded Summaries
#'
#' Builds a regular lattice over a point-occurrence data set, assigns every record
#' a `grid_id`, and returns grid-level summaries in both abundance and
#' presence/absence formats.
#'
#' The function:
#'   1. Tiles the study extent with square cells of size `grid_size`.
#'   2. Computes cell centroids and, for geographic CRS data, optional mapsheet
#'      codes (1:250 000 "BB" series).
#'   3. Aggregates user-specified columns (`sum_cols`) per cell, producing
#'      * **`grid_spp`**    - counts / abundances
#'      * **`grid_spp_pa`** - binary 0 / 1 table for dissimilarity analyses.
#'   4. Calculates helper metrics (`obs_sum`, `spp_rich`) for each cell.
#'   5. Rasterises key layers (`grid_id`, `obs_sum`, `spp_rich`) and preserves
#'      any extra metadata supplied via `extra_cols`.
#'
#' @param data            Data frame of points with x-y coordinates.
#' @param x_col,y_col    Column names for longitude and latitude.
#' @param grid_size       Cell size (degrees or metres, depending on CRS).
#' @param sum_cols        Character or numeric vector of columns to aggregate. Note: Numeric indices are converted to names internally.
#' @param extra_cols      Additional columns to keep (optional).
#' @param crs_epsg        EPSG code of the coordinate reference system.
#' @param unit            One of "deg", "min", "sec", or "m".
#'
#' @return A list containing
#'   - `grid_r`      - multi-layer **SpatRaster**
#'   - `grid_sf`     - polygon lattice with centroids & metrics
#'   - `grid_spp`    - abundance summary (*data.frame*)
#'   - `grid_spp_pa` - presence/absence summary (*data.frame*)
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' data = data.frame(
#'   x = runif(100, -10, 10),
#'   y = runif(100, -10, 10),
#'   species1 = rpois(100, 5),
#'   species2 = rpois(100, 3),
#'   recordedBy = sample(LETTERS, 100, replace = TRUE)
#' )
#' grid_result = generate_grid(data, x_col = "x", y_col = "y",
#'                              grid_size = 1, sum_cols = 3:4,
#'                              extra_cols = c("recordedBy"))
#' print(grid_result$block_sp)
#' plot(grid_result$grid_sf["grid_id"])
generate_grid <- function(data,
                          x_col      = "x",
                          y_col      = "y",
                          grid_size  = 0.5,
                          sum_cols   = NULL,      # char or numeric
                          extra_cols = NULL,
                          crs_epsg   = 4326,
                          unit       = c("deg", "min", "sec", "m"))
{
  # checks & setup
  unit <- match.arg(unit)
  stopifnot(all(c(x_col, y_col) %in% names(data)))

  if (is.null(sum_cols))
    stop("Please supply `sum_cols` (column names or indices of species).")

  # convert numeric indices -> names, verify all present
  if (is.numeric(sum_cols)) {
    if (max(sum_cols) > ncol(data))
      stop("`sum_cols` indices exceed number of columns in `data`.")
    sum_cols <- names(data)[sum_cols]
  }
  if (!all(sum_cols %in% names(data)))
    stop("Some columns in `sum_cols` not found in `data`.")

  sp_cols <- sum_cols                                    # use names henceforth

  # packages
  for (pkg in c("sf", "terra", "dplyr", "tidyr"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required but not installed.")

  # helper: presence/absence converter
  make_pa <- function(df)
    dplyr::mutate(df,
                  dplyr::across(dplyr::all_of(sp_cols), ~ ifelse(!is.na(.) & . > 0, 1, 0)))

  # save original coords & convert to sf
  data$orig_x <- data[[x_col]]
  data$orig_y <- data[[y_col]]
  pts_sf      <- sf::st_as_sf(data, coords = c(x_col, y_col), crs = crs_epsg)

  # branch 1 : grid_size == 0  ->  one site per unique coordinate
  if (grid_size == 0) {
    message("grid_size = 0 -> no lattice; grouping by unique coordinates.")

    # unique coords -> grid_id
    data <- dplyr::left_join(
      data,
      data.frame(orig_x = data$orig_x, orig_y = data$orig_y) |>
        dplyr::distinct() |>
        dplyr::mutate(grid_id = as.character(dplyr::row_number())),
      by = c("orig_x", "orig_y"))

    # summarise per grid_id
    grid_spp <- data |>
      dplyr::group_by(grid_id, dplyr::across(dplyr::all_of(extra_cols))) |>
      dplyr::summarise(dplyr::across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)),
                       .groups = "drop") |>
      dplyr::mutate(
        obs_sum  = rowSums(dplyr::pick(all_of(sp_cols)), na.rm = TRUE),
        spp_rich = rowSums(dplyr::pick(all_of(sp_cols)) > 0, na.rm = TRUE))

    grid_spp <- dplyr::relocate(
      grid_spp, grid_id, orig_x, orig_y, obs_sum, spp_rich)

    return(list(
      grid_r      = NULL,
      grid_sf     = NULL,
      grid_spp    = grid_spp,
      grid_spp_pa = make_pa(grid_spp)))
  }

  # branch 2 : build lattice
  message("Generating ", grid_size, "-", unit, " grid ...")

  # bounding box (expand a little)
  bb <- sf::st_bbox(pts_sf)
  bb[c("xmin","ymin")] <- floor(bb[c("xmin","ymin")] / grid_size) * grid_size -
    2 * grid_size
  bb[c("xmax","ymax")] <- ceiling(bb[c("xmax","ymax")] / grid_size) * grid_size +
    2 * grid_size

  # lattice polygons
  # Create grid polygons and compute centroids safely
  grid_polygons <- sf::st_make_grid(
    sf::st_as_sfc(bb), cellsize = grid_size, what = "polygons")

  grid_sf <- sf::st_sf(geometry = grid_polygons, crs = crs_epsg)

  # add centroid geometry and coordinates
  grid_sf$centroid     <- sf::st_centroid(grid_sf$geometry)
  coords               <- sf::st_coordinates(grid_sf$centroid)
  grid_sf$centroid_lon <- coords[, 1]
  grid_sf$centroid_lat <- coords[, 2]

  # supply unique grid IDs
  grid_sf$grid_id <- as.character(seq_len(nrow(grid_sf)))

  # optional 1:250k mapsheet code (deg/min only)
  if (unit %in% c("deg","min")) {
    lon_int <- floor(grid_sf$centroid_lon)
    lat_int <- floor(grid_sf$centroid_lat)
    grid_sf$mapsheet <- sprintf("%s%03d%s%02dBB",
                                ifelse(lon_int>=0,"E","W"), abs(lon_int),
                                ifelse(lat_int>=0,"N","S"), abs(lat_int))
  } else {
    grid_sf$mapsheet <- NA
  }

  # link points to cells
  data$grid_id <- sf::st_join(pts_sf, grid_sf["grid_id"],
                              join = sf::st_within)$grid_id

  # summarise species per cell
  grid_spp <- data |>
    dplyr::group_by(grid_id, dplyr::across(dplyr::all_of(extra_cols))) |>
    dplyr::summarise(dplyr::across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::mutate(
      obs_sum  = rowSums(dplyr::pick(all_of(sp_cols)), na.rm = TRUE),
      spp_rich = rowSums(dplyr::pick(all_of(sp_cols)) > 0, na.rm = TRUE))

  # attach centroids & mapsheets
  grid_spp <- dplyr::left_join(
    grid_spp,
    grid_sf |>
      sf::st_drop_geometry() |>
      dplyr::select(grid_id, centroid_lon, centroid_lat, mapsheet),
    by = "grid_id") |>
    dplyr::relocate(grid_id, centroid_lon, centroid_lat, mapsheet,
                    obs_sum, spp_rich)

  # add obs_sum & spp_rich back to polygons
  grid_sf <- dplyr::left_join(
    grid_sf, grid_spp |> dplyr::select(grid_id, obs_sum, spp_rich),
    by = "grid_id")

  # rasters: grid_id, obs_sum, spp_rich
  template <- terra::rast(xmin = bb["xmin"], xmax = bb["xmax"],
                          ymin = bb["ymin"], ymax = bb["ymax"],
                          resolution = grid_size,
                          crs = paste0("EPSG:", crs_epsg))
  grid_vect <- terra::vect(grid_sf)
  layers <- c("grid_id","obs_sum","spp_rich")
  rast_final <- do.call(c, lapply(layers, function(v){
    r <- terra::rasterize(grid_vect, template, field = v)
    names(r) <- v
    r
  }))

  # return
  list(
    grid_r      = rast_final,
    grid_sf     = grid_sf,
    grid_spp    = grid_spp,
    grid_spp_pa = make_pa(grid_spp))
}
