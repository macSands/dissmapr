#' Format Data Frame
#'
#' This function formats biodiversity data into long or wide format for
#' further analysis.
#'
#' @param data A `data.frame` of biodiversity records.
#' @param format Character; one of `"long"` or `"wide"`. If `NULL`, inferred
#'   from presence of `species_col` & `value_col`.
#' @param x_col Name of longitude column. If `NULL`, will search common names.
#' @param y_col Name of latitude column. If `NULL`, will search common names.
#' @param site_id_col Name of site identifier column. If `NULL`, will be
#'   generated from `(x, y)`.
#' @param species_col Name of species column (required for `"long"`).
#' @param value_col Name of value column (e.g. presence/abundance; for `"long"`).
#' @param sp_col_range Integer vector indexing species columns (for `"wide"`).
#' @param extra_cols Character vector of other columns to keep (e.g. metadata).
#'
#' @return A `list` with elements:
#'   - `site_obs`: data frame in **long** format (only for `format="long"`).
#'   - `site_sp`: data frame in **wide** format (species as columns).
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
    site_sp <- site_obs %>%
      dplyr::group_by(site_id, x, y, dplyr::across(dplyr::any_of(extra_cols)), species) %>%
      dplyr::summarize(value = sum(value, na.rm=TRUE), .groups="drop") %>%
      tidyr::pivot_wider(
        names_from  = species,
        values_from = value,
        values_fill = 0
      )

    # Make unique row names
    row.names(site_sp) <- make.unique(as.character(site_sp$site_id))

    return(list(site_obs = site_obs, site_sp = site_sp))
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

    site_sp <- data2 %>%
      dplyr::group_by(site_id, x, y, dplyr::across(dplyr::any_of(extra_cols))) %>%
      dplyr::summarize(
        dplyr::across(all_of(sp_cols), ~ sum(.x, na.rm=TRUE)),
        .groups = "drop"
      ) %>%
      as.data.frame()

    row.names(site_sp) <- make.unique(as.character(site_sp$site_id))

    return(list(site_sp = site_sp))
  }

  stop("`format` must be either 'long' or 'wide'")
}
