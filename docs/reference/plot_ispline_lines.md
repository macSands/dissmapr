# Plot ispline partial effects with quantile and start-point markers

Given a tidy data frame of ispline basis outputs (from
[`run_ispline_models`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md)),
this function identifies the ispline column matching a specified
covariate, draws thinner ispline curves for each `zOrder`, and overlays:

- Small symbols at user-defined quantiles of the raw covariate.

- A larger symbol at each curve's starting point (minimum covariate
  value).

## Usage

``` r
plot_ispline_lines(
  ispline_data,
  x_var,
  orders = NULL,
  cols = NULL,
  shapes = NULL,
  probs = c(0, 0.25, 0.5, 0.75, 1),
  line_size = 0.5,
  point_size = 1.5,
  start_size = 3,
  start_stroke = 0
)
```

## Arguments

- ispline_data:

  A data frame as returned by
  [`run_ispline_models`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md),
  containing raw covariates, their `_is` spline bases, and a `zOrder`
  column.

- x_var:

  A string; the name of the raw covariate to plot (e.g. “dist”). The
  function will look for a matching spline column named `<x_var>_is`.

- orders:

  Character vector of `zOrder` levels (in desired legend/order). Default
  is `unique(ispline_data$zOrder)`.

- cols:

  Character vector of colours, one per order. Default uses
  [`scales::hue_pal()`](https://scales.r-lib.org/reference/pal_hue.html).

- shapes:

  Integer vector of plotting symbols (pch codes), one per order. Default
  is `15:`(...).

- probs:

  Numeric vector of probabilities (between 0 and 1) at which to place
  small quantile markers. Default `c(0, .25, .5, .75, 1)`.

- line_size:

  Numeric; line width for the spline curves. Default 0.5.

- point_size:

  Numeric; size of the quantile markers. Default 1.5.

- start_size:

  Numeric; size of the big start-point markers. Default 3.

- start_stroke:

  Numeric; stroke width for the big start markers. Default 0.

## Value

A ggplot2 object showing:

- `geom_line()` for each curve.

- `geom_point()` at the specified quantiles.

- A larger `geom_point()` at each curve's minimum.

## Details

1.  Detects all columns ending in `_is` and fuzzy-matches `x_var`
    against them.

2.  Extracts the raw covariate name by stripping `_is`.

3.  Computes one "start" point per `zOrder` at the minimum raw
    covariate.

4.  Computes quantile-closest points per `zOrder` at the user's `probs`.

## See also

[`run_ispline_models`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# load sample data
data(bird.spec.fine)
data(bird.env.fine)

# prepare inputs
xy.bird    <- bird.spec.fine[,1:2]
spp.bird   <- bird.spec.fine[,3:102]
envir.bird <- bird.env.fine[,3:9]

# fit & gather ispline tables
ispline_tabs_all <- run_ispline_models(
  spp_df    = spp.bird,
  env_df    = envir.bird,
  xy_df     = xy.bird,
  orders    = 2:6,
  sam       = 100,
  normalize = "Jaccard",
  reg_type  = "ispline"
)

# Line plot for "dist"
plot_ispline_lines(
  ispline_data = ispline_tabs_all,
  x_var        = "dist",            # will match "dist_is"
  orders       = paste("Order", 2:6),
  cols         = c('green','cyan','purple','blue','black'),
  shapes       = c(15,16,17,18,19)
)
} # }
```
