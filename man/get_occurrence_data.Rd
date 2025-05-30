% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_occurrence_data.R
\name{get_occurrence_data}
\alias{get_occurrence_data}
\title{Import and harmonise biodiversity-occurrence data}
\usage{
get_occurrence_data(
  data = NULL,
  source_type = c("local_csv", "data_frame", "gbif"),
  gbif_zip_url = NULL,
  download_dir = tempdir(),
  sep = ",",
  site_id_col = NULL,
  x_col = NULL,
  y_col = NULL,
  sp_name_col = NULL,
  pa_col = NULL,
  abund_col = NULL,
  species_cols = NULL
)
}
\arguments{
\item{data}{File path (when \code{source_type = "local_csv"}),
an in-memory \code{data.frame} (\code{"data_frame"}), or \code{NULL}
(\code{"gbif"}).}

\item{source_type}{\code{"local_csv"}, \code{"data_frame"}, or \code{"gbif"}.}

\item{gbif_zip_url}{URL to a GBIF download ZIP (required when
\code{source_type = "gbif"}).}

\item{download_dir}{Folder to save the ZIP/extracted file
(default: \code{\link[=tempdir]{tempdir()}}).}

\item{sep}{Field separator for CSVs (default \code{","}).}

\item{site_id_col, x_col, y_col, sp_name_col, pa_col, abund_col}{Optional custom column names.}

\item{species_cols}{\strong{Optional} numeric or character vector specifying
the species columns in a wide input (e.g. \code{4:11} or
\code{c("Sp1","Sp2")}).  Overrides the default \code{"sp_*"}
detection.}
}
\value{
A \code{data.frame}:
\describe{
\item{Long format}{Columns \code{site_id}, \code{x}, \code{y}, \code{sp_name}, plus \code{pa}
\emph{or} \code{abund}.}
\item{Wide → long}{Same columns after stacking the specified or
auto-detected species columns.}
}
}
\description{
\code{get_occurrence_data()} reads occurrence records from a \strong{local CSV/file
path}, an \strong{in-memory \code{data.frame}}, or a \strong{GBIF download (ZIP)} and
returns a tidy data frame in either \strong{long} (one row = one record) or
\strong{wide} (one row = one site, one column = one species) form.
}
\details{
Column names are auto-detected from common patterns
(\code{"site_id"}, \code{"x"}, \code{"y"}, \code{"sp_name"}, \code{"pa"} or \code{"abund"}).
Supply \verb{*_col} arguments \strong{only} when your data use different names.

For wide data the helper normally looks for columns that start with
\code{"sp_"}.  Set \code{species_cols} to a numeric range (e.g. \code{4:11}) or a character
vector of column names when the species columns do \strong{not} follow the
\code{"sp_*"} convention.
}
\section{Workflow}{

\enumerate{
\item \strong{Read} the data from \code{source_type}.
\item \strong{Detect / insert} compulsory columns (site, coords, species, value).
\item \strong{Validate} coordinates (-180 ≤ lon ≤ 180, -90 ≤ lat ≤ 90).
\item \strong{Return}
\itemize{
\item a long table (\code{site_id}, \code{x}, \code{y}, \code{sp_name}, \code{pa|abund}) when species
name + value columns are present; or
\item a long table reshaped from wide species columns.
}
}
}

\examples{
# 1. Local CSV example -----------------------------------------------
tmp <- tempfile(fileext = ".csv")
df_local <- data.frame(
  site_id = 1:10,
  x = runif(10), y = runif(10),
  sp_name = c("plant1", "plant2","plant3", "plant4","plant5", "plant1", "plant2","plant3", "plant4","plant5"),
  abun = sample(0:20, size = 10, replace = TRUE)
)
write.csv(df_local, tmp, row.names = FALSE)
local_test = get_occurrence_data(data = tmp, source_type = "local_csv", sep = ",")

# 2. Existing wide-format data.frame -----------------------------------------------
df_wide <- df_local \%>\%
  pivot_wider(
    names_from  = sp_name,   # these become column names
    values_from = abun,  # fill these cell values
    values_fn   = sum,       # sum duplicates
    values_fill = 0          # fill missing with 0
  )
wide_test = get_occurrence_data(data = df_wide, source_type = "data_frame", species_cols = 4:11)

# 3. Custom names ----------------------------------------------------------
names(sim_dat)[1:5] <- c("plot_id", "lon", "lat", "taxon", "presence")
occ_long2 <- get_occurrence_data(
  data           = sim_dat,
  source_type    = "data_frame",
  site_id_col    = "plot_id",
  x_col          = "lon",
  y_col          = "lat",
  sp_name_col    = "taxon",
  pa_col         = "presence"
)
head(occ_long2)

# 4. GBIF download (requires internet) -----------------------------------------------
\dontrun{
gbif_test = get_occurrence_data(
  source_type   = "gbif",
  gbif_zip_url  = "https://api.gbif.org/v1/occurrence/download/request/0038969-240906103802322.zip"
)
}
}
\seealso{
\code{\link[tidyr:pivot_longer]{tidyr::pivot_longer()}} used internally.
}
