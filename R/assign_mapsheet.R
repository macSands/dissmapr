# #==============================================================================
# Author: Sandra MacFadyen
# Email: sandra@biogis.co.za
# Date: 20 Dec 2024
#==============================================================================

# More info on mapsheet codes:
# https://ngi.dalrrd.gov.za/index.php/what-we-do/maps-and-geospatial-information/41-sa-mapsheet-referencing
# https://www.schoolgrades.co.za/mapwork-grade-12-notes-geography-study-guides/


#' Add Nearest Mapsheet Code and Center Coordinates
#'
#' This function takes an existing data frame with coordinate columns and
#' appends:
#'  1. the nearest mapsheet code (format "[E|W]DDD[N|S]DDBB"), and
#'  2. the center longitude/latitude of that mapsheet.
#'
#' It divides the world into regular grid cells (in degrees, minutes, seconds,
#' or metres), finds which cell each point falls into, and then formats the
#' sheet code based on the cell centre.
#'
#' @param data A data.frame containing at least the coordinate columns.
#' @param x_col Character. Name of the longitude (x) column. Default `"x"`.
#' @param y_col Character. Name of the latitude (y) column. Default `"y"`.
#' @param cell_size Numeric. Cell size: e.g. 1 for 1° cells, 15 for 15′ cells.
#' @param unit Character. Unit of `cell_size`; one of:
#'   - `"deg"` for decimal degrees (default)
#'   - `"min"` for arc-minutes
#'   - `"sec"` for arc-seconds
#'   - `"m"`   for metres (projected coords).
#'
#' @return The input `data` with three new columns:
#'   - `mapsheet`: the 7-character sheet code
#'   - `center_lon`: centre longitude of that sheet
#'   - `center_lat`: centre latitude of that sheet
#'
#' @examples
#' bird_obs = data.frame(
#'   site_id = 1:6,
#'   x = c(22.71862, 20.40034, 18.51368, 18.38477, 23.56160, 18.87285),
#'   y = c(-33.98912, -34.45408, -34.08271, -34.25946, -33.97620, -34.06225),
#'   species = c("Fulica cristata", "Numida meleagris coronatus",
#'               "Anas undulata", "Oenanthe familiaris",
#'               "Cyanomitra veroxii", "Gallinula chloropus"),
#'   value = rep(1, 6),
#'   year = c(2023, 2022, 2023, 2016, 2016, 2019),
#'   month = c(5, 12, 10, 8, 9, 2),
#'   day = c(12, 4, 3, 1, 4, 21)
#' )
#'
#' # For 1° mapsheet cells:
#' bird_obs2 = add_nearest_mapsheet(bird_obs, cell_size = 1, unit = "deg")
#' head(bird_obs2)
#'
#' @export
assign_mapsheet = function(data,
                            x_col = "x",
                            y_col = "y",
                            cell_size = 1,
                            unit = c("deg", "min", "sec", "m")) {
  unit = match.arg(unit)

  # Check that coordinate columns exist
  if (!x_col %in% names(data)) {
    stop("Column '", x_col, "' not found in the data.")
  }
  if (!y_col %in% names(data)) {
    stop("Column '", y_col, "' not found in the data.")
  }

  # Convert cell_size to coordinate units (degrees or equivalent)
  if (unit == "deg") {
    inc = cell_size
  } else if (unit == "min") {
    inc = cell_size / 60
  } else if (unit == "sec") {
    inc = cell_size / 3600
  } else if (unit == "m") {
    inc = cell_size  # For projected coordinates (e.g. meters)
  }

  # Compute the center of the grid cell in which each point falls.
  # The lower boundary is given by floor(coord / inc)*inc,
  # so the cell center is lower + (inc/2).
  data$center_lon = floor(data[[x_col]] / inc) * inc + inc/2
  data$center_lat = floor(data[[y_col]] / inc) * inc + inc/2

  # Generate mapsheet code using the center coordinates.
  # The mapsheet code format is:
  #   [E|W] + zero-padded 3-digit integer part of center_lon +
  #   [N|S] + zero-padded 2-digit integer part of center_lat +
  #   sub-cell codes (here fixed as "B" for both).
  lon_int = floor(data$center_lon)
  lat_int = floor(data$center_lat)
  lon_dir = ifelse(data$center_lon >= 0, "E", "W")
  lat_dir = ifelse(data$center_lat >= 0, "N", "S")
  data$mapsheet = sprintf("%s%03d%s%02d%s%s",
                           lon_dir, abs(lon_int),
                           lat_dir, abs(lat_int),
                           "B", "B")

  return(data)
}
