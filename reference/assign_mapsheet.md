# Add Nearest Mapsheet Code and Center Coordinates

This function takes an existing data frame with coordinate columns and
appends:

1.  the nearest mapsheet code (format `"[E|W]DDD[N|S]DDBB"`), and

2.  the center longitude/latitude of that mapsheet.

## Usage

``` r
assign_mapsheet(
  data,
  x_col = "x",
  y_col = "y",
  cell_size = 1,
  unit = c("deg", "min", "sec", "m")
)
```

## Arguments

- data:

  A data.frame containing at least the coordinate columns.

- x_col:

  Character. Name of the longitude (x) column. Default `"x"`.

- y_col:

  Character. Name of the latitude (y) column. Default `"y"`.

- cell_size:

  Numeric. Cell size: e.g. 1 for 1deg cells, 15 for 15' cells.

- unit:

  Character. Unit of `cell_size`; one of:

  - `"deg"` for decimal degrees (default)

  - `"min"` for arc-minutes

  - `"sec"` for arc-seconds

  - `"m"` for metres (projected coords).

## Value

The input `data` with three new columns:

- `mapsheet`: the 7-character sheet code

- `center_lon`: centre longitude of that sheet

- `center_lat`: centre latitude of that sheet

## Details

It divides the world into regular grid cells (in degrees, minutes,
seconds, or metres), finds which cell each point falls into, and then
formats the sheet code based on the cell centre.

## Examples

``` r
bird_obs = data.frame(
  site_id = 1:6,
  x = c(22.71862, 20.40034, 18.51368, 18.38477, 23.56160, 18.87285),
  y = c(-33.98912, -34.45408, -34.08271, -34.25946, -33.97620, -34.06225),
  species = c("Fulica cristata", "Numida meleagris coronatus",
              "Anas undulata", "Oenanthe familiaris",
              "Cyanomitra veroxii", "Gallinula chloropus"),
  value = rep(1, 6),
  year = c(2023, 2022, 2023, 2016, 2016, 2019),
  month = c(5, 12, 10, 8, 9, 2),
  day = c(12, 4, 3, 1, 4, 21)
)

# For 1deg mapsheet cells:
bird_obs2 = assign_mapsheet(bird_obs, cell_size = 1, unit = "deg")
head(bird_obs2)
#>   site_id        x         y                    species value year month day
#> 1       1 22.71862 -33.98912            Fulica cristata     1 2023     5  12
#> 2       2 20.40034 -34.45408 Numida meleagris coronatus     1 2022    12   4
#> 3       3 18.51368 -34.08271              Anas undulata     1 2023    10   3
#> 4       4 18.38477 -34.25946        Oenanthe familiaris     1 2016     8   1
#> 5       5 23.56160 -33.97620         Cyanomitra veroxii     1 2016     9   4
#> 6       6 18.87285 -34.06225        Gallinula chloropus     1 2019     2  21
#>   center_lon center_lat  mapsheet
#> 1       22.5      -33.5 E022S34BB
#> 2       20.5      -34.5 E020S35BB
#> 3       18.5      -34.5 E018S35BB
#> 4       18.5      -34.5 E018S35BB
#> 5       23.5      -33.5 E023S34BB
#> 6       18.5      -34.5 E018S35BB
```
