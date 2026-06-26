# Calculate Distance Between Sites

Calculate Distance Between Sites

## Usage

``` r
distance(df, site_col, vec_from, vec_to = NULL, coord_cols = c("x", "y"))
```

## Arguments

- df:

  A data frame containing site coordinates.

- site_col:

  The column name representing site IDs.

- vec_from:

  The site ID or vector for the starting site.

- vec_to:

  The site ID or vector for the destination site(s).

## Value

Distance(s) in meters between the specified sites.
