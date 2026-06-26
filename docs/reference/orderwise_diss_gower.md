# Compute Gower dissimilarity between two site vectors

Calculates the Gower dissimilarity (a bounded 0-1 measure) between two
numeric vectors (e.g., species abundances) using the
[`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html)
function. When the second vector is `NULL` (order 1), returns `NA_real_`
since a pairwise comparison is not meaningful.

## Usage

``` r
orderwise_diss_gower(vec_from, vec_to = NULL)
```

## Arguments

- vec_from:

  Numeric vector representing the first site's attributes (e.g., species
  counts or trait values).

- vec_to:

  Optional numeric vector for the second site; if `NULL`, the function
  returns `NA_real_`.

## Value

A single numeric dissimilarity value between 0 and 1, or `NA_real_` if
`vec_to` is `NULL`.
