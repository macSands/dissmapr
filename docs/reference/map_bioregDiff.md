# Map Bioregional Change Metrics Between Categorical Raster Layers

Calculates five complementary indices that quantify how the categorical
(e.g. bioregion / cluster) label of each raster cell changes across a
temporal or scenario stack of rasters:

- **Difference count** - total number of times a cell's label differs
  from the first layer.

- **Shannon entropy** - information-theoretic diversity of labels within
  the cell's time-series.

- **Stability** - proportion of layers in which the label is identical
  to the first layer *(1 = always unchanged, 0 = always different)*.

- **Transition frequency** - sum of binary change maps between
  successive layers (*how often a change occurs between any pair of
  neighbours*).

- **Weighted change index** - cumulative dissimilarity-weighted change
  where the weight is derived from the empirical frequency of
  transitions between all pairs of labels.

## Usage

``` r
map_bioregDiff(raster_input, approach = "all")
```

## Arguments

- raster_input:

  A multi-layer `SpatRaster` **or** a `list` of single-layer
  `SpatRaster` objects representing the same spatial extent/resolution.

- approach:

  Character string specifying the metric to return:
  `"difference_count"`, `"shannon_entropy"`, `"stability"`,
  `"transition_frequency"`, `"weighted_change_index"`, or `"all"`
  (default) for a five-layer stack containing every metric.

## Value

A `SpatRaster`:

- **single-layer** if `approach` is one of the named metrics;

- **five-layer** (names: `Difference_Count`, `Shannon_Entropy`,
  `Stability`, `Transition_Frequency`, `Weighted_Change_Index`) if
  `approach = "all"`.

## Details

The dissimilarity weights for the **weighted change index** are built
from the observed transition table of successive layers, normalised to
lie between 0 and 1 (larger = rarer transition). The function accepts
either a multi-layer `SpatRaster` or a plain `list` of single-layer
`SpatRaster`s, which is internally concatenated with **terra**.

## See also

[`app`](https://rspatial.github.io/terra/reference/app.html),
[`rast`](https://rspatial.github.io/terra/reference/rast.html)

## Examples

``` r
## -------------------------------------------------------------
## Minimal reproducible example with three random categorical
## rasters (four classes, identical geometry)
## -------------------------------------------------------------
if (requireNamespace("terra", quietly = TRUE)) {
  set.seed(42)

  r1 <- terra::rast(nrows = 40, ncols = 40,
                    vals  = sample(1:4, 40 * 40, TRUE))
  r2 <- terra::rast(r1, vals = sample(1:4, terra::ncell(r1), TRUE))
  r3 <- terra::rast(r1, vals = sample(1:4, terra::ncell(r1), TRUE))

  r_stack <- terra::rast(list(r1, r2, r3))
  names(r_stack) <- paste0("t", 1:3)

  ## 1. All five metrics
  diff_all <- map_bioregDiff(r_stack, approach = "all")
  print(diff_all)

  ## 2. Just the Shannon-entropy layer
  ent <- map_bioregDiff(r_stack, approach = "shannon_entropy")
  terra::plot(ent, main = "Shannon entropy of label sequence")
}
#> class       : SpatRaster
#> size        : 40, 40, 5  (nrow, ncol, nlyr)
#> resolution  : 9, 4.5  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84)
#> source(s)   : memory
#> names       : Differ~_Count, Shanno~ntropy, Stability, Transi~quency, Weight~_Index
#> min values  :             0,            -0,         0,             0,             0
#> max values  :             2,      1.098612,         1,             2,      0.325991

```
