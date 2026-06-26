# Run multiple Zeta.msgdm ispline models and return both models and combined ispline table

Fits `Zeta.msgdm` models of type “ispline” for a series of zeta-orders,
extracts the raw environmental covariates (plus distance) and their
ispline bases, and returns both the list of fitted models and one tidy
data frame combining all orders.

## Usage

``` r
run_ispline_models(
  spp_df,
  env_df,
  xy_df,
  orders = 2:6,
  sam = 100,
  distance.type = "Euclidean",
  normalize = "Jaccard",
  reg_type = "ispline"
)
```

## Arguments

- spp_df:

  A data frame or matrix of species incidence/abundance.

- env_df:

  A data frame of environmental covariates.

- xy_df:

  A two-column data frame or matrix of site coordinates.

- orders:

  Integer vector of zeta orders to fit (e.g. 2:6).

- sam:

  Integer; number of random samples per order (passed to `Zeta.msgdm`).

- distance.type:

  Character; distance metric for `Zeta.msgdm` (default “Euclidean”).

- normalize:

  Character; normalization method for `Zeta.msgdm` (default “Jaccard”).

- reg_type:

  Character; regression type for `Zeta.msgdm` (default “ispline”).

## Value

A named list with:

- `zeta_gdm_list`:

  A list of the fitted `Zeta.msgdm()` objects, named by order.

- `ispline_table`:

  A tibble with one row per sample, containing all original covariates
  (including `distance`), the ispline bases suffixed `_is`, and a
  `zOrder` column.

## See also

[`Zeta.msgdm`](https://rdrr.io/pkg/zetadiv/man/Zeta.msgdm.html),
[`Return.ispline`](https://rdrr.io/pkg/zetadiv/man/Return.ispline.html)

## Examples

``` r
if (FALSE) { # \dontrun{
data(bird.spec.fine); data(bird.env.fine)
xy   <- bird.spec.fine[,1:2]
spp  <- bird.spec.fine[,3:102]
env  <- bird.env.fine[,3:9]

out <- run_ispline_models(
  spp_df        = spp,
  env_df        = env,
  xy_df         = xy,
  orders        = 2:6,
  sam           = 100,
  normalize     = "Jaccard",
  reg_type      = "ispline"
)
names(out)
head(out$ispline_table)
} # }
```
