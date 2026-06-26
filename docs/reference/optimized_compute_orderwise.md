# Compute Order-wise Pairwise Metrics with Optional Parallel Sampling

A flexible, high-performance engine for calculating dissimilarity,
distance, or any custom metric between a *focal* site and every
combination of \\(order-1)\\ other sites. The routine supports:

- special fast paths for `order = 1` (trivial) and `order = 2` when
  `func = orderwise_diss_gower`;

- great-circle distances via `geodist_helper` (requires `coord_cols`);

- arbitrary metrics acting on numeric "site vectors" defined by
  `sp_cols`;

- parallel execution with **future.apply** and progress bars from
  **pbapply**;

- random or complete sampling of higher-order combinations.

## Usage

``` r
optimized_compute_orderwise(
  df,
  func,
  site_col,
  sp_cols,
  order = 2,
  sample_no = NULL,
  sample_portion = 1,
  parallel = TRUE,
  n_workers = parallel::detectCores() - 1,
  coord_cols = NULL
)
```

## Arguments

- df:

  Data frame or data.table containing the site data.

- func:

  Function of two numeric vectors returning a single numeric value (e.g.
  `orderwise_diss_gower`, `geodist_helper`).

- site_col:

  Character; column with the unique site identifiers.

- sp_cols:

  Character or integer vector giving the columns used to build each
  site's numeric vector **when `func` is not `geodist_helper`**.

- order:

  Integer scalar or vector of orders (\>= 1). For values \\\> 2\\ the
  function iterates through each order.

- sample_no:

  Integer. Maximum *number* of combinations to sample **per focal site**
  for orders \>= 3. If `NULL` (default) all combinations are used
  (subject to `sample_portion`).

- sample_portion:

  Numeric in (0, 1\]. Proportion of combinations to retain when
  `sample_no = NULL`. Ignored when `sample_no` is supplied.

- parallel:

  Logical. Run computations in parallel with **future.apply** (default
  `TRUE`).

- n_workers:

  Integer. Number of background workers; default
  `parallel::detectCores() - 1`.

- coord_cols:

  Character vector of coordinate columns (e.g. `c("x","y")`) required
  when `func = geodist_helper`.

## Value

A `data.table` with four columns

- site_from:

  Focal site ID.

- site_to:

  Comma-separated list of comparison sites (`NA` when `order = 1`).

- value:

  Computed metric.

- order:

  The order used for the calculation.

## Details

For each site the relevant data vector is pre-extracted once and
re-used, avoiding repeated sub-setting. When `order = 2` and
`func = orderwise_diss_gower` the full pairwise Gower distance matrix is
built in one call to
[`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html),
melted, and returned, yielding a major speed-up over looping.

## Dependencies

Requires **pbapply**, **data.table**, **future.apply**, **reshape2**,
and **dplyr** (checked at run-time). For Gower distances the example
assumes
[`orderwise_diss_gower()`](https://b-cubed-eu.github.io/dissmapr/reference/orderwise_diss_gower.md)
is exported by *dissmapr*.

## Examples

``` r
## ---------------------------------------------------------------
## Example 1 - Gower dissimilarity between species vectors (order 2)
## ---------------------------------------------------------------
if (requireNamespace("cluster", quietly = TRUE)) {
  set.seed(1)
  n_sites  <- 8
  n_spp    <- 20
  test_df  <- data.frame(
    site_id = paste0("s", seq_len(n_sites)),
    x       = runif(n_sites, 23, 24),
    y       = runif(n_sites, -34, -33),
    matrix(rpois(n_sites * n_spp, 1),
           nrow = n_sites,
           dimnames = list(NULL, paste0("sp", seq_len(n_spp))))
  )

  res2 <- optimized_compute_orderwise(
    df        = test_df,
    func      = orderwise_diss_gower,   # dissmapr helper
    site_col  = "site_id",
    sp_cols   = paste0("sp", 1:n_spp),
    order     = 2,
    parallel  = FALSE
  )
  head(res2)
}
#> Error in optimized_compute_orderwise(df = test_df, func = orderwise_diss_gower,     site_col = "site_id", sp_cols = paste0("sp", 1:n_spp), order = 2,     parallel = FALSE): could not find function "optimized_compute_orderwise"

## ---------------------------------------------------------------
## Example 2 - Great-circle distance (km) between site centroids
## ---------------------------------------------------------------
if (requireNamespace("geosphere", quietly = TRUE)) {
  # geodist_helper is expected to be exported by dissmapr
  res_geo <- optimized_compute_orderwise(
    df         = test_df,
    func       = geodist_helper,
    site_col   = "site_id",
    order      = 2,
    coord_cols = c("x", "y"),
    parallel   = FALSE
  )
  head(res_geo)
}
#> Error in optimized_compute_orderwise(df = test_df, func = geodist_helper,     site_col = "site_id", order = 2, coord_cols = c("x", "y"),     parallel = FALSE): could not find function "optimized_compute_orderwise"
```
