# Plot facetted boxplots for all ispline basis columns

Creates a multi-facet boxplot of every ispline basis column in your
data, grouped by a specified order factor (e.g. zeta orders). Each
spline term (columns ending in `_is`) gets its own facet, and uses a
color-blind-friendly Viridis palette.

## Usage

``` r
plot_ispline_boxplots(
  ispline_data,
  ispline_suffix = "_is",
  order_col = "zOrder",
  palette = "viridis",
  direction = -1,
  ncol = 3,
  outlier_size = 0.5
)
```

## Arguments

- ispline_data:

  A data frame as returned by
  [`run_ispline_models`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md),
  containing raw covariates, their spline bases (suffixed `_is`), and an
  order column.

- ispline_suffix:

  A string suffix identifying your spline columns. Default is “\_is”.

- order_col:

  The name of the grouping column (e.g. “zOrder”).

- palette:

  One of the Viridis options (“viridis”, “magma”, “plasma”, “cividis”,
  etc.). Default “viridis”.

- direction:

  Integer 1 or -1 to control palette direction. Default -1 (reversed).

- ncol:

  Number of columns in the facet wrap. Default 3.

- outlier_size:

  Size of the outlier points. Default 0.5.

## Value

A ggplot2 object showing one boxplot per spline term.

## Details

1.  Automatically detects all columns whose names end with
    `ispline_suffix`.

2.  Pivots the data to long format for ggplot2.

3.  Facets a boxplot for each spline variable with independent scales.

## See also

[`run_ispline_models`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md),
[`plot_ispline_lines`](https://b-cubed-eu.github.io/dissmapr/reference/plot_ispline_lines.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# load example data
data(bird.spec.fine)
data(bird.env.fine)

# prepare inputs
xy.bird    <- bird.spec.fine[,1:2]
spp.bird   <- bird.spec.fine[,3:102]
envir.bird <- bird.env.fine[,3:9]

# Fit & gather ispline tables
ispline_tabs_all <- run_ispline_models(
  spp_df    = spp.bird,
  env_df    = envir.bird,
  xy_df     = xy.bird,
  orders    = 2:6,
  sam       = 100,
  normalize = "Jaccard",
  reg_type  = "ispline"
)

# Facetted boxplots of all *_is columns
plot_ispline_boxplots(
  ispline_data   = ispline_tabs_all,
  ispline_suffix = "_is",
  order_col      = "zOrder",
  palette        = "viridis",
  direction      = -1,
  ncol           = 3
)
} # }
```
