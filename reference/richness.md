# Calculate Species Richness

Counts unique species per site; reports absolute differences for
pairwise and summary stats (range, variance, mean) for multi-site.

## Usage

``` r
richness(vec_from, vec_to = NULL)
```

## Arguments

- vec_from:

  A numeric vector representing species counts at the first site.

- vec_to:

  (Optional) A numeric vector for pairwise or higher-order comparisons.

## Value

The species richness value.
