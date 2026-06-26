# Format Biodiversity Records to Long / Wide

Converts a table of biodiversity observations between **long** and
**wide** layouts while standardising key column names.

- **Long format** - one row per observation with columns `site_id`, `x`,
  `y`, `species`, `value` (+ optional `extra_cols`).

- **Wide format** - one row per site with species as individual columns.

## Usage

``` r
format_df(
  data,
  format = NULL,
  x_col = NULL,
  y_col = NULL,
  site_id_col = NULL,
  species_col = NULL,
  value_col = NULL,
  sp_col_range = NULL,
  extra_cols = NULL
)
```

## Arguments

- data:

  A data frame containing biodiversity records.

- format:

  Character; target layout `"long"` or `"wide"`. If `NULL` the format is
  inferred automatically.

- x_col, y_col:

  Character. Names of the longitude (*x*) and latitude (*y*) columns. If
  `NULL`, common alternatives are searched.

- site_id_col:

  Character. Column giving a unique site identifier. If `NULL`, a new
  `site_id` is generated from the coordinate pair.

- species_col:

  Character. Column containing species names (required for
  `format = "long"`).

- value_col:

  Character. Column with numeric values such as presence/absence (0/1)
  or abundance. If `NULL`, each record is assigned a value of 1.

- sp_col_range:

  Integer vector giving the index of species columns when
  `format = "wide"`. If `NULL` all non-coordinate / non-metadata columns
  are treated as species.

- extra_cols:

  Character vector of additional columns to carry through to the output
  (e.g. sampling metadata or environmental covariates).

## Value

A named list with up to two elements

- `site_obs` - a long-format data frame (returned only when
  `format = "long"`).

- `site_spp` - a wide site x species data frame.

## Details

If column names are not supplied, the function attempts to detect common
variants (e.g. `"lon"`, `"longitude"` for *x*). When converting long -\>
wide, duplicate observations of the same species at a site are
aggregated by summing `value`. When converting wide -\> long, species
columns are inferred either from `sp_col_range` or by excluding
coordinate / metadata columns.

## Dependencies

Relies on **dplyr**, **tidyr**, and **rlang** (loaded with
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)).

## See also

[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html),
[`pivot_wider`](https://tidyr.tidyverse.org/reference/pivot_wider.html)

## Examples

``` r
## --- Example 1: long  ->  wide --------------------------------------------
ex_long <- data.frame(
  lon     = c(23.10, 23.10, 23.25, 23.25),
  lat     = c(-34.00, -34.00, -34.05, -34.05),
  species = c("sp1",  "sp2",  "sp1",  "sp3"),
  count   = c(1, 2, 3, 1)
)

out_long <- format_df(
  data        = ex_long,
  format      = "long",
  x_col       = "lon",
  y_col       = "lat",
  species_col = "species",
  value_col   = "count"
)

head(out_long$site_spp)
#>        site_id     x      y sp1 sp2 sp3
#> site_1  site_1 23.10 -34.00   1   2   0
#> site_2  site_2 23.25 -34.05   3   0   1

## --- Example 2: wide ->  long --------------------------------------------
ex_wide <- out_long$site_spp

out_wide <- format_df(
  data   = ex_wide,
  format = "wide"
)

head(out_wide$site_spp)
#>        site_id     x      y sp1 sp2 sp3
#> site_1  site_1 23.10 -34.00   1   2   0
#> site_2  site_2 23.25 -34.05   3   0   1
```
