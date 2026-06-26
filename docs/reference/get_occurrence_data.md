# Import and harmonise biodiversity-occurrence data

`get_occurrence_data()` reads occurrence records from a **local CSV/file
path**, an **in-memory `data.frame`**, or a **GBIF download (ZIP)** and
returns a tidy data frame in either **long** (one row = one record) or
**wide** (one row = one site, one column = one species) form.

## Usage

``` r
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
```

## Arguments

- data:

  File path (when `source_type = "local_csv"`), an in-memory
  `data.frame` (`"data_frame"`), or `NULL` (`"gbif"`).

- source_type:

  `"local_csv"`, `"data_frame"`, or `"gbif"`.

- gbif_zip_url:

  URL to a GBIF download ZIP (required when `source_type = "gbif"`).

- download_dir:

  Folder to save the ZIP/extracted file (default:
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- sep:

  Field separator for CSVs (default `","`).

- site_id_col, x_col, y_col, sp_name_col, pa_col, abund_col:

  Optional custom column names.

- species_cols:

  **Optional** numeric or character vector specifying the species
  columns in a wide input (e.g. `4:11` or `c("Sp1","Sp2")`). Overrides
  the default `"sp_*"` detection.

## Value

A `data.frame`:

- Long format:

  Columns `site_id`, `x`, `y`, `sp_name`, plus `pa` *or* `abund`.

- Wide - long:

  Same columns after stacking the specified or auto-detected species
  columns.

## Details

Column names are auto-detected from common patterns (`"site_id"`, `"x"`,
`"y"`, `"sp_name"`, `"pa"` or `"abund"`). Supply `*_col` arguments
**only** when your data use different names.

For wide data the helper normally looks for columns that start with
`"sp_"`. Set `species_cols` to a numeric range (e.g. `4:11`) or a
character vector of column names when the species columns do **not**
follow the `"sp_*"` convention.

## Workflow

1.  **Read** the data from `source_type`.

2.  **Detect / insert** compulsory columns (site, coords, species,
    value).

3.  **Validate** coordinates (-180 \<= lon \<= 180, -90 \<= lat \<= 90).

4.  **Return**

    - a long table (`site_id`, `x`, `y`, `sp_name`, `pa|abund`) when
      species name + value columns are present; or

    - a long table reshaped from wide species columns.

## See also

[`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)
used internally.

## Examples

``` r
# 1. Long-format data read from a local CSV file ---------------------------
tmp <- tempfile(fileext = ".csv")
df_local <- data.frame(
  site_id = rep(1:5, each = 2),
  x = runif(10, 18, 20), y = runif(10, -34, -33),
  sp_name = rep(c("sp_a", "sp_b"), times = 5),
  abund = sample(0:20, size = 10, replace = TRUE)
)
write.csv(df_local, tmp, row.names = FALSE)
local_test <- get_occurrence_data(data = tmp, source_type = "local_csv")
head(local_test)
#>   site_id        x         y sp_name abund
#> 1       1 18.07474 -33.18352    sp_a    11
#> 2       1 19.03761 -33.96050    sp_b    18
#> 3       2 19.35803 -33.26097    sp_a     3
#> 4       2 19.80647 -33.65128    sp_b     2
#> 5       3 18.05105 -33.17075    sp_a    13
#> 6       3 19.97816 -33.46446    sp_b     7

# 2. Wide-format data.frame (species as columns) ---------------------------
df_wide <- data.frame(
  site_id = 1:5,
  x = runif(5, 18, 20), y = runif(5, -34, -33),
  sp_a = sample(0:5, 5, replace = TRUE),
  sp_b = sample(0:5, 5, replace = TRUE),
  sp_c = sample(0:5, 5, replace = TRUE)
)
wide_test <- get_occurrence_data(
  data = df_wide, source_type = "data_frame", species_cols = 4:6
)

# 3. Custom column names ---------------------------------------------------
sim_dat <- data.frame(
  plot_id  = 1:6,
  lon      = runif(6, 18, 20),
  lat      = runif(6, -34, -33),
  taxon    = rep(c("sp_a", "sp_b"), 3),
  presence = 1
)
occ_long2 <- get_occurrence_data(
  data        = sim_dat,
  source_type = "data_frame",
  site_id_col = "plot_id",
  x_col       = "lon",
  y_col       = "lat",
  sp_name_col = "taxon",
  pa_col      = "presence"
)
head(occ_long2)
#>   site_id        x         y sp_name presence pa
#> 1       1 18.06394 -33.19857    sp_a        1  1
#> 2       2 18.47446 -33.85372    sp_b        1  1
#> 3       3 19.37298 -33.17728    sp_a        1  1
#> 4       4 18.45164 -33.66900    sp_b        1  1
#> 5       5 18.63699 -33.62583    sp_a        1  1
#> 6       6 18.34797 -33.37025    sp_b        1  1

# 4. GBIF download (requires internet) -------------------------------------
if (FALSE) { # \dontrun{
gbif_test <- get_occurrence_data(
  source_type  = "gbif",
  gbif_zip_url = "https://api.gbif.org/v1/occurrence/download/request/0038969-240906103802322.zip"
)
} # }
```
