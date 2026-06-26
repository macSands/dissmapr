# Predict Pairwise Compositional Turnover (zeta-dissimilarity) with Richness

Takes raw species and environmental data, fits a multi-site GDM model
output, computes predicted pairwise turnover (zeta2) across the
landscape, and returns a data frame with site-level richness,
environmental covariates, distance, and predicted turnover; optionally
plots a heatmap of zeta2 predictions.

## Usage

``` r
predict_dissim(
  grid_spp,
  species_cols,
  env_vars,
  zeta_model,
  grid_xy,
  bndy_fc = NULL,
  x_col = "x",
  y_col = "y",
  show_plot = TRUE,
  skip_scale = FALSE
)
```

## Arguments

- grid_spp:

  Data frame containing site IDs, coordinates, and species
  presence-absence/abundance columns.

- species_cols:

  Integer or character vector giving the columns of `grid_spp` that hold
  species data.

- env_vars:

  Data frame of raw environmental predictors (unscaled; rows must align
  with `grid_spp`).

- zeta_model:

  Fitted object from
  [`Zeta.msgdm`](https://rdrr.io/pkg/zetadiv/man/Zeta.msgdm.html) (order
  = 2, reg.type = "ispline").

- grid_xy:

  Data frame of site coordinates, same row-order as `grid_spp`, with
  columns `x_col`, `y_col`.

- bndy_fc:

  Optional `sf` or `SpatVector` polygon to overlay.

- x_col:

  Name of the x (longitude) column in `grid_spp`/`grid_xy`.

- y_col:

  Name of the y (latitude) column.

- show_plot:

  Logical; if TRUE (default), print the turnover heatmap.

- skip_scale:

  Logical; if TRUE, skip internal scaling/standardisation of the
  environmental predictors (use when `env_vars` are already scaled).
  Default FALSE.

## Value

A data frame with one row per site, containing:

- richness:

  Species richness (sum across `species_cols`).

- distance:

  Mean great-circle distance (km) from each site to all others.

- \<env_vars\>:

  All scaled environmental predictors.

- pred_zeta:

  Linear predictor (logit scale) from `Predict.msgdm()`.

- pred_zetaExp:

  Predicted turnover (0-1 scale).

- log_pred_zetaExp:

  Natural log of `pred_zetaExp`.

- x_col, y_col:

  Site coordinates (from `grid_xy`).

## Examples

``` r
if (FALSE) { # \dontrun{
result <- predict_dissim(
  grid_spp     = bird.spec.fine,
  species_cols = 3:102,
  env_vars     = bird.env.fine[,3:9],
  zeta_model   = z_mod,
  grid_xy      = bird.spec.fine[,1:2],
  bndy_fc      = rsa,
  x_col        = "x",
  y_col        = "y",
  show_plot    = FALSE
)
} # }
```
