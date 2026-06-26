# Compute Order-wise Metrics

This function computes metrics for ecological data across specified
order levels. It supports single-site, pairwise, and higher-order
calculations, and allows for parallel processing for efficiency.

## Usage

``` r
compute_orderwise(
  df,
  func,
  site_col,
  sp_cols = NULL,
  order = 2,
  sample_no = NULL,
  sample_portion = 1,
  parallel = TRUE,
  n_workers = parallel::detectCores() - 1
)
```

## Arguments

- df:

  A data frame containing the ecological data.

- func:

  A function to compute metrics. It must accept inputs in the form of
  species vectors or site information depending on the order.

- site_col:

  A character string specifying the column name in `df` representing
  site IDs.

- sp_cols:

  A vector of column names in `df` representing species data (default:
  NULL).

- order:

  An integer or vector of integers specifying the order(s) of
  computation.

  - `1`: Single-site computations.

  - `2`: Pairwise computations.

  - `>= 3`: Higher-order computations.

- sample_no:

  An integer specifying the maximum number of combinations to sample for
  higher-order computations (default: NULL for all combinations).

- sample_portion:

  A numeric value between 0 and 1 indicating the proportion of
  combinations to sample for higher-order computations (default: 1,
  meaning 100%).

- parallel:

  A logical value indicating whether to enable parallel computation
  (default: TRUE).

- n_workers:

  An integer specifying the number of parallel workers to use (default:
  one less than the number of available cores).

## Value

A `data.table` containing the results of computations. Columns include:

- `site_from`: The source site.

- `site_to`: The target site(s) (NA for order = 1).

- `order`: The computation order.

- `value`: The computed metric value.

## Examples

``` r
# Example usage with a custom metric function
metric_func = function(vec1, vec2 = NULL) {
  if (is.null(vec2)) {
    return(sum(vec1))  # Example: species richness
  } else {
    return(sum((vec1 - vec2)^2))  # Example: pairwise dissimilarity
  }
}

data = data.frame(
  site = rep(letters[1:3], each = 3),
  sp1 = c(1, 0, 2, 1, 2, 0, 0, 1, 1),
  sp2 = c(0, 1, 1, 2, 0, 0, 1, 0, 2)
)

# SPECIES RICHNESS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
richness = function(vec_from, vec_to = NULL) {
 if (is.null(vec_to)) {
   # Handle single-site calculations (order = 1)
   return(sum(vec_from != 0, na.rm = TRUE))
 } else if (length(vec_from) > 1 && length(vec_to) > 1) {
   # Handle pairwise or higher-order comparisons
   return(abs(sum(vec_from != 0, na.rm = TRUE) - sum(vec_to != 0, na.rm = TRUE)))
 } else {
   # Invalid input case
   return(NA)
 }
}

rich_o12 = compute_orderwise(
 df = data,
 func = richness,
 site_col = 'site',
 sp_cols = c('sp1', 'sp2'),
 order = 1:2,      # Compute for pairwise and higher-order comparisons
 parallel = FALSE  # keep the example single-threaded
)
#> Time elapsed for order 1: 0 minutes and 0.00 seconds
#> Time elapsed for order 2: 0 minutes and 0.00 seconds
#> Total computation time: 0 minutes and 0.00 seconds
head(rich_o12)
#>    site_from site_to order value
#>       <char>  <char> <int> <int>
#> 1:         a    <NA>     1    11
#> 2:         b    <NA>     1    11
#> 3:         c    <NA>     1    11
#> 4:         b       a     2     0
#> 5:         c       a     2     0
#> 6:         a       b     2     0
```
