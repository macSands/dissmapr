#' Format Biodiversity Records to Long / Wide
#'
#' @description
#' Converts a table of biodiversity observations between **long** and
#' **wide** layouts while standardising key column names.
#' * **Long format** – one row per observation with columns
#'   `site_id`, `x`, `y`, `species`, `value` (+ optional `extra_cols`).
#' * **Wide format** – one row per site with species as individual columns.
#'
#' @details
#' If column names are not supplied, the function attempts to detect common
#' variants (e.g. `"lon"`, `"longitude"` for *x*).
#' When converting long → wide, duplicate observations of the same species at a
#' site are aggregated by summing `value`.
#' When converting wide → long, species columns are inferred either from
#' `sp_col_range` or by excluding coordinate / metadata columns.
#'
#' @param data A data frame containing biodiversity records.
#' @param format Character; target layout `"long"` or `"wide"`.
#'   If `NULL` the format is inferred automatically.
#' @param x_col,y_col Character. Names of the longitude (*x*) and latitude (*y*)
#'   columns. If `NULL`, common alternatives are searched.
#' @param site_id_col Character. Column giving a unique site identifier.
#'   If `NULL`, a new `site_id` is generated from the coordinate pair.
#' @param species_col Character. Column containing species names
#'   (required for `format = "long"`).
#' @param value_col Character. Column with numeric values such as
#'   presence/absence (0/1) or abundance. If `NULL`, each record is assigned
#'   a value of 1.
#' @param sp_col_range Integer vector giving the index of species columns when
#'   `format = "wide"`. If `NULL` all non-coordinate / non-metadata columns are
#'   treated as species.
#' @param extra_cols Character vector of additional columns to carry through to
#'   the output (e.g. sampling metadata or environmental covariates).
#'
#' @return A named list with up to two elements
#' * `site_obs` – a long-format data frame (returned only when
#'   `format = "long"`).
#' * `site_spp` – a wide site × species data frame.
#'
#' @section Dependencies:
#' Relies on **dplyr**, **tidyr**, and **rlang** (loaded with
#' `requireNamespace()`).
#'
#' @seealso \code{\link[dplyr]{group_by}}, \code{\link[tidyr]{pivot_wider}}
#'
#' @importFrom dplyr group_by across mutate rename select summarise filter any_of
#'   cur_group_id ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom rlang sym
#'
#' @examples
#' ## --- Example 1: long  →  wide --------------------------------------------
#' ex_long <- data.frame(
#'   lon     = c(23.10, 23.10, 23.25, 23.25),
#'   lat     = c(-34.00, -34.00, -34.05, -34.05),
#'   species = c("sp1",  "sp2",  "sp1",  "sp3"),
#'   count   = c(1, 2, 3, 1)
#' )
#'
#' out_long <- format_df(
#'   data        = ex_long,
#'   format      = "long",
#'   x_col       = "lon",
#'   y_col       = "lat",
#'   species_col = "species",
#'   value_col   = "count"
#' )
#'
#' head(out_long$site_spp)
#'
#' ## --- Example 2: wide →  long --------------------------------------------
#' ex_wide <- out_long$site_spp
#'
#' out_wide <- format_df(
#'   data   = ex_wide,
#'   format = "wide"
#' )
#'
#' head(out_wide$site_spp)
#'
#' @export
format_df <- function(data,
                      format       = NULL,
                      x_col        = NULL,
                      y_col        = NULL,
                      site_id_col  = NULL,
                      species_col  = NULL,
                      value_col    = NULL,
                      sp_col_range = NULL,
                      extra_cols   = NULL) {
  # Dependencies
  if (!requireNamespace("dplyr", quietly=TRUE) ||
      !requireNamespace("tidyr", quietly=TRUE) ||
      !requireNamespace("rlang", quietly=TRUE)) {
    stop("Please install dplyr, tidyr, and rlang first.")
  }

  # Null-coalesce operator
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Helper to find a column from alternatives
  find_col <- function(alts) {
    nm <- names(data)
    m  <- tolower(nm) %in% tolower(alts)
    if (any(m)) nm[which(m)[1]] else NULL
  }

  # Resolve column names (or detect defaults)
  x_col       <- x_col       %||% find_col(c("x","lon","longitude","decimalLongitude"))
  y_col       <- y_col       %||% find_col(c("y","lat","latitude","decimalLatitude"))
  site_id_col <- site_id_col %||% find_col(c("site_id","grid_id","id","plot"))
  species_col <- species_col %||% find_col(c("species","sp_name","verbatimScientificName"))
  value_col   <- value_col   %||% find_col(c("pa","presence","abundance","count"))

  stopifnot(!is.null(x_col), !is.null(y_col))

  # Create site_id if missing
  if (is.null(site_id_col)) {
    data <- data %>%
      dplyr::group_by(dplyr::across(all_of(c(x_col,y_col)))) %>%
      dplyr::mutate(site_id = paste0("site_", dplyr::cur_group_id())) %>%
      dplyr::ungroup()
    site_id_col <- "site_id"
  }

  # Drop any existing 'species' column if renaming another column to 'species'
  if (!is.null(species_col) &&
      !identical(species_col, "species") &&
      "species" %in% names(data)) {
    data <- data[, setdiff(names(data), "species"), drop = FALSE]
  }

  # Infer format if not provided
  format <- format %||%
    if (!is.null(species_col) && !is.null(value_col)) "long" else "wide"

  ## ---- LONG format -----------------------------------
  if (format == "long") {
    stopifnot(!is.null(species_col))

    # Standardize names via tidy-eval
    data2 <- data %>%
      dplyr::rename(
        site_id = !!rlang::sym(site_id_col),
        x       = !!rlang::sym(x_col),
        y       = !!rlang::sym(y_col),
        species = !!rlang::sym(species_col)
      ) %>%
      dplyr::mutate(
        value = if (!is.null(value_col))
          as.numeric(!!rlang::sym(value_col))
        else 1
      ) %>%
      dplyr::filter(!is.na(species) & species != "")

    # Build site_obs
    site_obs <- data2 %>%
      dplyr::select(site_id, x, y, species, value, dplyr::any_of(extra_cols))

    # Pivot to wide (one column per species)
    site_spp <- site_obs %>%
      dplyr::group_by(site_id, x, y, dplyr::across(dplyr::any_of(extra_cols)), species) %>%
      dplyr::summarize(value = sum(value, na.rm=TRUE), .groups="drop") %>%
      tidyr::pivot_wider(
        names_from  = species,
        values_from = value,
        values_fill = 0
      )

    # Make unique row names
    row.names(site_spp) <- make.unique(as.character(site_spp$site_id))

    return(list(site_obs = site_obs, site_spp = site_spp))
  }

  ## ---- WIDE format -----------------------------------
  if (format == "wide") {
    # Determine species columns
    if (!is.null(sp_col_range)) {
      sp_cols <- names(data)[sp_col_range]
    } else {
      sp_cols <- setdiff(names(data),
                         c(site_id_col, x_col, y_col, extra_cols))
    }
    stopifnot(length(sp_cols) > 0)

    data2 <- data %>%
      dplyr::rename(
        site_id = !!rlang::sym(site_id_col),
        x       = !!rlang::sym(x_col),
        y       = !!rlang::sym(y_col)
      )

    site_spp <- data2 %>%
      dplyr::group_by(site_id, x, y, dplyr::across(dplyr::any_of(extra_cols))) %>%
      dplyr::summarize(
        dplyr::across(all_of(sp_cols), ~ sum(.x, na.rm=TRUE)),
        .groups = "drop"
      ) %>%
      as.data.frame()

    row.names(site_spp) <- make.unique(as.character(site_spp$site_id))

    return(list(site_spp = site_spp))
  }

  stop("`format` must be either 'long' or 'wide'")
}
