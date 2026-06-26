# Retrieve, crop, resample, and link environmental rasters to sampling sites

`get_enviro_data()` automates six steps:

1.  **AOI build** - buffers the convex hull of all input points.

2.  **Raster acquisition** - downloads or loads a multi-layer stack
    (WorldClim, SoilGrids, footprint, population, or user-supplied).

3.  **Cropping** - trims the stack to the buffered AOI.

4.  **Optional resampling** - resamples the cropped stack to a
    user-defined grid resolution.

5.  **Extraction & gap fill** - pulls raster values at each site and
    linearly interpolates isolated NAs.

6.  **Assembly** - returns a tidy *site x environment* table plus the
    cropped (and resampled) raster and an `sf` layer of sites.

## Usage

``` r
get_enviro_data(
  data,
  buffer_km = 10,
  source = "geodata",
  var = "bio",
  res = 2.5,
  grid_r = NULL,
  path = "data/",
  year = NULL,
  depth = NULL,
  stat = "mean",
  model = NULL,
  ssp = NULL,
  time = NULL,
  sp_cols = NULL,
  ext_cols = NULL
)
```

## Arguments

- data:

  Data frame of spatial points; must include lon/lat columns such as
  `"x","y"` or `"centroid_lon","centroid_lat"`.

- buffer_km:

  Buffer width (km) for the AOI. Default 10.

- source:

  `"geodata"` (default) to fetch layers via **geodata** or `"local"` to
  supply a local `SpatRaster` or file path.

- var:

  Raster product to download (see details) or ignored when
  `source = "local"`.

- res:

  Resolution (arc-min) for WorldClim/WorldPop layers (geodata).

- grid_r:

  Optional grid resolution (in same CRS units) to which the cropped
  raster is resampled before extraction. If `NULL`, no resampling is
  performed.

- path:

  Cache folder for downloaded rasters (created if absent).

- year, model, ssp, time:

  Optional arguments for time-stamped products (human footprint,
  population, CMIP6, etc.).

- depth, stat:

  Arguments passed to
  [`geodata::soil_world()`](https://rspatial.github.io/geodata/reference/soil_grids.html).

- sp_cols:

  **Columns to drop** from the final table (e.g. a large species
  matrix). Accepts names or numeric indices *relative to `data`*.

- ext_cols:

  **Columns to append** verbatim (e.g. `"obs_sum","spp_rich"`).

## Value

A list with

- `env_rast` `SpatRaster` - cropped (and optionally resampled)
  environmental stack

- `sites_sf` `sf` POINT layer (WGS-84) of the input sites

- `env_df` Tibble with site ID, coordinates, every raster variable, plus
  any columns requested in `ext_cols`

## Examples

``` r
if (FALSE) { # \dontrun{
# Example sites with longitude/latitude columns
sites <- data.frame(
  site_id = 1:3,
  x = c(18.5, 19.0, 19.5),
  y = c(-33.9, -33.5, -34.1)
)

# Download WorldClim bioclimatic variables and link them to the sites
enviro <- get_enviro_data(
  data   = sites,
  source = "geodata",
  var    = "bio",
  res    = 5
)
head(enviro$env_df)
} # }
```
