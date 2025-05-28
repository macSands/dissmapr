#' Import and harmonise biodiversity-occurrence data
#'
#' `get_occurrence_data()` reads occurrence records from a **local CSV/file
#' path**, an **in-memory `data.frame`**, or a **GBIF download (ZIP)** and
#' returns a tidy data frame in either **long** (one row = one record) or
#' **wide** (one row = one site, one column = one species) form.
#'
#' Column names are auto-detected from common patterns
#' (`"site_id"`, `"x"`, `"y"`, `"sp_name"`, `"pa"` or `"abund"`).
#' Supply `*_col` arguments **only** when your data use different names.
#'
#' For wide data the helper normally looks for columns that start with
#' `"sp_"`.  Set `species_cols` to a numeric range (e.g. `4:11`) or a character
#' vector of column names when the species columns do **not** follow the
#' `"sp_*"` convention.
#'
#' @section Workflow:
#' 1. **Read** the data from `source_type`.
#' 2. **Detect / insert** compulsory columns (site, coords, species, value).
#' 3. **Validate** coordinates (-180 ≤ lon ≤ 180, -90 ≤ lat ≤ 90).
#' 4. **Return**
#'    * a long table (`site_id`, `x`, `y`, `sp_name`, `pa|abund`) when species
#'      name + value columns are present; or
#'    * a long table reshaped from wide species columns.
#'
#' @param data           File path (when `source_type = "local_csv"`),
#'                       an in-memory `data.frame` (`"data_frame"`), or `NULL`
#'                       (`"gbif"`).
#' @param source_type    `"local_csv"`, `"data_frame"`, or `"gbif"`.
#' @param gbif_zip_url   URL to a GBIF download ZIP (required when
#'                       `source_type = "gbif"`).
#' @param download_dir   Folder to save the ZIP/extracted file
#'                       (default: [tempdir()]).
#' @param sep            Field separator for CSVs (default `","`).
#' @param site_id_col,x_col,y_col,sp_name_col,pa_col,abund_col
#'                       Optional custom column names.
#' @param species_cols   **Optional** numeric or character vector specifying
#'                       the species columns in a wide input (e.g. `4:11` or
#'                       `c("Sp1","Sp2")`).  Overrides the default `"sp_*"`
#'                       detection.
#'
#' @return A `data.frame`:
#' \describe{
#'   \item{Long format}{Columns `site_id`, `x`, `y`, `sp_name`, plus `pa`
#'     *or* `abund`.}
#'   \item{Wide → long}{Same columns after stacking the specified or
#'     auto-detected species columns.}
#' }
#'
#' @seealso [tidyr::pivot_longer()] used internally.
#'
#' @importFrom dplyr rename mutate left_join
#' @importFrom tidyr pivot_longer
#' @importFrom httr GET write_disk
#' @importFrom data.table fread
#' @importFrom utils unzip
#' @export
#'
#' @examples
#' # 1. Local CSV example -----------------------------------------------
#' tmp <- tempfile(fileext = ".csv")
#' df_local <- data.frame(
#'   site_id = 1:10,
#'   x = runif(10), y = runif(10),
#'   sp_name = c("plant1", "plant2","plant3", "plant4","plant5", "plant1", "plant2","plant3", "plant4","plant5"),
#'   abun = sample(0:20, size = 10, replace = TRUE)
#' )
#' write.csv(df_local, tmp, row.names = FALSE)
#' local_test = get_occurrence_data(data = tmp, source_type = "local_csv", sep = ",")
#'
#' # 2. Existing wide-format data.frame -----------------------------------------------
#' df_wide <- df_local %>%
#'   pivot_wider(
#'     names_from  = sp_name,   # these become column names
#'     values_from = abun,  # fill these cell values
#'     values_fn   = sum,       # sum duplicates
#'     values_fill = 0          # fill missing with 0
#'   )
#' wide_test = get_occurrence_data(data = df_wide, source_type = "data_frame", species_cols = 4:11)
#'
#' # 3. Custom names ----------------------------------------------------------
#' names(sim_dat)[1:5] <- c("plot_id", "lon", "lat", "taxon", "presence")
#' occ_long2 <- get_occurrence_data(
#'   data           = sim_dat,
#'   source_type    = "data_frame",
#'   site_id_col    = "plot_id",
#'   x_col          = "lon",
#'   y_col          = "lat",
#'   sp_name_col    = "taxon",
#'   pa_col         = "presence"
#' )
#' head(occ_long2)
#'
#' # 4. GBIF download (requires internet) -----------------------------------------------
#' \dontrun{
#' gbif_test = get_occurrence_data(
#'   source_type   = "gbif",
#'   gbif_zip_url  = "https://api.gbif.org/v1/occurrence/download/request/0038969-240906103802322.zip"
#' )
#' }
get_occurrence_data <- function(
    data           = NULL,
    source_type    = c("local_csv", "data_frame", "gbif"),
    gbif_zip_url   = NULL,
    download_dir   = tempdir(),
    sep            = ",",
    site_id_col    = NULL,
    x_col          = NULL,
    y_col          = NULL,
    sp_name_col    = NULL,
    pa_col         = NULL,
    abund_col      = NULL,
    species_cols   = NULL
) {

  ## ---- 0.  dependencies ------------------------------------------------ ##
  pkgs <- c("dplyr", "tidyr", "httr", "data.table")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
  }

  source_type <- match.arg(source_type)

  ## helper: detect a column name (case-insensitive) ---------------------- ##
  detect_columns <- function(dat, candidates) {
    hits <- intersect(tolower(names(dat)), tolower(candidates))
    if (length(hits))
      return(names(dat)[tolower(names(dat)) %in% hits][1L])
    NULL
  }

  ## ---- 1.  read data --------------------------------------------------- ##
  if (source_type == "local_csv") {
    if (!is.character(data) || !file.exists(data))
      stop("File not found: '", data, "'.")
    data <- utils::read.csv(
      data, sep = sep, stringsAsFactors = FALSE,
      row.names = NULL, check.names = FALSE
    )

  } else if (source_type == "data_frame") {
    if (!is.data.frame(data))
      stop("'data' must be a data.frame when source_type = 'data_frame'.")

  } else {                              # ---- GBIF ----------------------- #
    if (is.null(gbif_zip_url))
      stop("'gbif_zip_url' is required for source_type = 'gbif'.")
    if (!dir.exists(download_dir))
      dir.create(download_dir, recursive = TRUE)

    zip_path <- file.path(download_dir, basename(gbif_zip_url))
    httr::GET(gbif_zip_url, httr::write_disk(zip_path, overwrite = TRUE))

    csv_list <- utils::unzip(zip_path, list = TRUE)
    occ_file <- csv_list$Name[
      grepl("(occurrence\\.txt|\\.csv)$", csv_list$Name, ignore.case = TRUE)
    ]
    if (!length(occ_file))
      stop("No occurrence file (*.csv or occurrence.txt) found in the ZIP.")

    utils::unzip(zip_path, files = occ_file, exdir = download_dir)
    data <- data.table::fread(
      file.path(download_dir, occ_file),
      sep = "\t", data.table = FALSE
    )
  }

  ## ---- 2.  candidate column names ------------------------------------- ##
  column_sets <- list(
    site_id = c("site_id", "site", "sample", "id", "plot"),
    x       = c("x", "lon", "long", "longitude", "x_coord", "decimalLongitude"),
    y       = c("y", "lat", "latitude", "y_coord", "decimalLatitude"),
    sp_name = c("sp_name", "name", "species", "scientific",
                "spp", "verbatimScientificName"),
    pa      = c("pa", "presence", "obs"),
    abund   = c("abund", "abundance", "count", "total")
  )

  ## prepend user names so they win detection ---------------------------- ##
  if (!is.null(site_id_col)) column_sets$site_id <- c(site_id_col, column_sets$site_id)
  if (!is.null(x_col))       column_sets$x       <- c(x_col,       column_sets$x)
  if (!is.null(y_col))       column_sets$y       <- c(y_col,       column_sets$y)
  if (!is.null(sp_name_col)) column_sets$sp_name <- c(sp_name_col, column_sets$sp_name)
  if (!is.null(pa_col))      column_sets$pa      <- c(pa_col,      column_sets$pa)
  if (!is.null(abund_col))   column_sets$abund   <- c(abund_col,   column_sets$abund)

  detected <- lapply(column_sets, detect_columns, dat = data)

  ## ---- 3.  coordinate validation -------------------------------------- ##
  if (is.null(detected$x) || is.null(detected$y))
    stop("Longitude and/or latitude columns not detected; ",
         "specify them with 'x_col' and 'y_col'.")

  if (any(data[[detected$x]] <  -180 | data[[detected$x]] >  180, na.rm = TRUE))
    stop("Longitude values must be between -180 and 180.")
  if (any(data[[detected$y]] <   -90 | data[[detected$y]] >   90, na.rm = TRUE))
    stop("Latitude values must be between -90 and 90.")

  ## ---- 4.  create site_id if missing ---------------------------------- ##
  if (is.null(detected$site_id)) {
    coords_unique <- unique(data[, c(detected$x, detected$y), drop = FALSE])
    coords_unique$site_id <- seq_len(nrow(coords_unique))
    data <- dplyr::left_join(
      data, coords_unique, by = c(detected$x, detected$y)
    )
    detected$site_id <- "site_id"
  }

  ## ---- 5.  ensure a value column -------------------------------------- ##
  if (is.null(detected$pa) && is.null(detected$abund)) {
    data <- dplyr::mutate(data, pa = 1)
    detected$pa <- "pa"
  }

  ## ---- 6.  decide format & reshape ------------------------------------ ##
  ## 6a. species columns in a wide table --------------------------------- ##
  if (!is.null(species_cols)) {
    if (is.numeric(species_cols)) {
      if (any(species_cols < 1 | species_cols > ncol(data)))
        stop("'species_cols' indices out of range.")
      sp_cols <- names(data)[species_cols]
    } else if (is.character(species_cols)) {
      if (!all(species_cols %in% names(data)))
        stop("Some names in 'species_cols' not found in data.")
      sp_cols <- species_cols
    } else {
      stop("'species_cols' must be numeric indices or character names.")
    }
  } else {
    sp_cols <- grep("^sp_", names(data), value = TRUE)
  }

  long_ok <- !is.null(detected$sp_name) &&
    (!is.null(detected$pa) || !is.null(detected$abund))
  wide_ok <- length(sp_cols) > 0L

  # ------------------ long --------- #
  if (long_ok) {
    data <- dplyr::rename(
      data,
      site_id = detected$site_id,
      x       = detected$x,
      y       = detected$y,
      sp_name = detected$sp_name
    )

    if (!is.null(detected$pa)) {
      data <- dplyr::mutate(
        data,
        pa = suppressWarnings(
          as.numeric(ifelse(is.na(.data[[detected$pa]]) |
                              .data[[detected$pa]] == "", 1, .data[[detected$pa]]))
        )
      )
    } else {
      data <- dplyr::mutate(
        data,
        abund = suppressWarnings(as.numeric(.data[[detected$abund]]))
      )
    }
    return(data)

    # ------------------ wide -------- #
  } else if (wide_ok) {
    data <- tidyr::pivot_longer(
      data,
      cols      = dplyr::all_of(sp_cols),  # <- wrap the vector
      names_to  = "sp_name",
      values_to = "value"
    )
    return(data)

  } else {
    stop("Unable to determine whether the input is long or wide; ",
         "provide species columns via 'species_cols' or ensure they start with 'sp_'.")
  }
}
