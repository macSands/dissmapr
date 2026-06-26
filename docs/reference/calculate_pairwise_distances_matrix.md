# Calculate Pairwise Distance Matrix

Calculate Pairwise Distance Matrix

## Usage

``` r
calculate_pairwise_distances_matrix(
  data,
  x_col = "x",
  y_col = "y",
  distance_fun = distGeo
)
```

## Arguments

- data:

  A data frame containing site coordinates and IDs.

- x_col:

  Name of the longitude column in `data` (default: "x").

- y_col:

  Name of the latitude column in `data` (default: "y").

- distance_fun:

  The distance function to use (default: `distGeo`).

## Value

A data frame containing pairwise distances between sites (in km).

## Examples

``` r
# default usage when your coords are in columns "x" and "y"
my_sites_df <- data.frame(
  grid_id = 1:4,
  x = c(18.4, 18.6, 19.0, 19.2),
  y = c(-33.9, -34.0, -33.7, -33.5)
)
dist_df_default <- calculate_pairwise_distances_matrix(data = my_sites_df)
head(dist_df_default)
#>   site_from site_to    value
#> 2         2       1 21.55999
#> 3         3       1 59.82529
#> 4         4       1 86.42350
#> 5         1       2 21.55999
#> 7         3       2 49.77614
#> 8         4       2 78.52510

# specify custom column names for longitude and latitude
sites_custom <- data.frame(
  grid_id = 1:4,
  centroid_lon = c(18.4, 18.6, 19.0, 19.2),
  centroid_lat = c(-33.9, -34.0, -33.7, -33.5)
)
dist_df <- calculate_pairwise_distances_matrix(
  data  = sites_custom,
  x_col = "centroid_lon",
  y_col = "centroid_lat"
)
head(dist_df)
#>   site_from site_to    value
#> 2         2       1 21.55999
#> 3         3       1 59.82529
#> 4         4       1 86.42350
#> 5         1       2 21.55999
#> 7         3       2 49.77614
#> 8         4       2 78.52510
```
