# Calculate geographic distance via Haversine formula (helper)

Computes the geographic (Haversine) distance in meters between two
coordinate pairs, or returns 0 for single-site (order=1) calculations
when `vec_to` is `NULL`.

## Usage

``` r
geodist_helper(
  vec_from,
  vec_to = NULL,
  coord_cols = c("centroid_lon", "centroid_lat")
)
```

## Arguments

- vec_from:

  Numeric vector of length 2, or a named vector containing coordinates.
  If named, should include names given by `coord_cols`.

- vec_to:

  Optional numeric vector of length 2 (destination coordinates); if
  `NULL` (default), the function returns 0.

- coord_cols:

  Character vector of length 2 giving the names of the longitude and
  latitude elements in `vec_from`/`vec_to` when those are named vectors.
  Defaults to `c("x", "y")`.

## Value

Numeric distance in meters between the two points, or 0 if `vec_to` is
`NULL`.
