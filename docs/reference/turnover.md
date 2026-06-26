# Calculate Species Turnover or Beta Diversity

Computes beta diversity as proportion unshared; 0-1 for pairwise,
extended to multi-site averages or totals.

## Usage

``` r
turnover(vec_from, vec_to = NULL)
```

## Arguments

- vec_from:

  A numeric vector representing species counts at the first site.

- vec_to:

  (Optional) A numeric vector for pairwise or higher-order comparisons.

## Value

The species turnover value.
