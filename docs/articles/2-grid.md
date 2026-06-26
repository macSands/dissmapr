# User-defined grid

## `dissmapr`

### A Novel Framework for Automated Compositional Dissimilarity and Biodiversity Turnover Analysis

#### 1. User-defined area of interest and grid resolution

Defining the geographic extent and an analysis grid early ensures that
all subsequent data extraction, aggregation, and visualisation tasks are
carried out within a consistent spatial framework. In this vignette we:

1.  **Load the national boundary of South Africa** to set our area of
    interest (AoI).
2.  **Select a working resolution** of 0.5° (≈ 55 km) to balance spatial
    detail with computational cost.
3.  **Convert the AoI to a `terra` vector** so that raster operations
    run efficiently.
4.  **Create a blank raster template** using the chosen resolution and
    the AoI’s CRS (Coordinate Reference System).
5.  **Populate the raster with placeholder values** (here simply 1).
6.  **Mask the raster to the AoI** so that only cells whose centroids
    fall within South Africa remain.

``` r
# 1. Load the national boundary 
# The shapefile is shipped with the package for full reproducibility.
rsa = sf::st_read(system.file("extdata", "rsa.shp", package = "dissmapr"))
#> Reading layer `rsa' from data source 
#>   `C:\Users\macfadyen\AppData\Local\R\cache\R\renv\library\dissmapr-ee60f76d\windows\R-4.5\x86_64-w64-mingw32\dissmapr\extdata\rsa.shp' 
#>   using driver `ESRI Shapefile'
#> Simple feature collection with 1 feature and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 16.45802 ymin: -34.83514 xmax: 32.89125 ymax: -22.12661
#> Geodetic CRS:  WGS 84

# 2. Choose a working resolution 
# A 0.5‑degree cell size strikes a balance between computational load and
# the spatial resolution at which national‑level biodiversity patterns remain
# interpretable.
res = 0.5   # decimal degrees° (≈ 55 km at the equator)

# 3. Convert the AoI to a 'terra' vector 
# 'terra' supports fast raster operations; converting now avoids repeated
# coercion later.
rsa_vect = terra::vect(rsa)

# 4. Initialise a blank raster template 
# The template inherits the AoI’s coordinate reference system (CRS) and is
# discretised into equally‑sized cells according to the resolution chosen.
grid = terra::rast(rsa_vect, resolution = res, crs = terra::crs(rsa_vect))

# 5. Populate the raster with placeholder values 
# We simply assign the value 1 to every cell; the values themselves are
# irrelevant at this stage—the grid’s geometry is what matters.
terra::values(grid) = 1

# 6. Clip the raster to the AoI 
# Any cells whose centroids fall outside the boundary are set to NA, thereby
# restricting subsequent computations to the AoI only.
grid_masked = terra::mask(grid, rsa_vect)
# grid_masked is now a 0.5° lattice clipped to South Africa and will serve as the common spatial denominator for all downstream summaries.
```

#### 2. Summarise records by grid centroid using `generate_grid()`

With the national lattice in place, we can now **condense point-level
observations to grid cells** using
[`generate_grid()`](https://b-cubed-eu.github.io/dissmapr/reference/generate_grid.md)
to:

1.  **Construct a bounding grid**: Expands the extent of input points
    and tessellates it with square cells of the chosen size (here 0.5°).
2.  **Allocate a `grid_id`**: Every record inherits the ID of the cell
    in which it falls.
3.  **Aggregate user-selected columns** within each occupied cell,
    returning:
    - `grid_spp`: species counts / abundances.  
    - `grid_spp_pa`: the same matrix recoded to presence (1) /
      absence (0) for binary dissimilarity metrics.  
    - `obs_sum`: total observations across the aggregated columns.  
    - `spp_rich`: number of columns with a non-zero count (simple
      species richness).
4.  **Compute cell centroids** and optional assign mapsheet codes
    (useful for atlasing projects).
5.  **Rasterise key layers** (`grid_id`, `obs_sum`, `spp_rich`) for fast
    map algebra.
6.  **Return four spatial objects** ready for further analysis:
    - `grid_r`: multi-layer `SpatRaster`  
    - `grid_sf`: polygon lattice with centroids & summaries  
    - `grid_spp`: abundance table (per cell × species)  
    - `grid_spp_pa`: binary presence/absence table (same dimensions as
      `grid_spp`)

Because every observation is now referenced to a regular grid, all
downstream statistics and graphics are standardised to the same sample
area.

``` r
# Generate a 0.5° grid summary for the point dataset `site_spp`
grid_list = dissmapr::generate_grid(
  data          = site_spp,           # point data with x/y + species columns
  x_col         = "x",                # longitude column
  y_col         = "y",                # latitude  column
  grid_size     = 0.5,                # cell size in degrees
  sum_cols      = 4:ncol(site_spp),   # columns to aggregate
  crs_epsg      = 4326                # WGS84
)

# Inspect the returned list 
str(grid_list, max.level = 1)
#> List of 4
#>  $ grid_r     :S4 class 'SpatRaster' [package "terra"]
#>  $ grid_sf    :Classes 'sf' and 'data.frame':    1110 obs. of  8 variables:
#>   ..- attr(*, "sf_column")= chr "geometry"
#>   ..- attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA NA NA NA NA
#>   .. ..- attr(*, "names")= chr [1:6] "centroid_lon" "centroid_lat" "grid_id" "mapsheet" ...
#>  $ grid_spp   : tibble [415 × 2,874] (S3: tbl_df/tbl/data.frame)
#>  $ grid_spp_pa: tibble [415 × 2,874] (S3: tbl_df/tbl/data.frame)

# (Optional) Promote list items to named objects 
grid_r = grid_list$grid_r   # raster summary
grid_sf = grid_list$grid_sf   # polygons for mapping or joins
grid_spp = grid_list$grid_spp # tabular summary per cell
grid_spp_pa = grid_list$grid_spp_pa # presence/absence summary

# Quick checks 
dim(grid_sf); head(grid_sf)
#> [1] 1110    8
#> Simple feature collection with 6 features and 6 fields
#> Active geometry column: geometry
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 15.5 ymin: -36 xmax: 18.5 ymax: -35.5
#> Geodetic CRS:  WGS 84
#>   centroid_lon centroid_lat grid_id  mapsheet obs_sum spp_rich
#> 1        15.75       -35.75       1 E015S36BB      NA       NA
#> 2        16.25       -35.75       2 E016S36BB      NA       NA
#> 3        16.75       -35.75       3 E016S36BB      NA       NA
#> 4        17.25       -35.75       4 E017S36BB      NA       NA
#> 5        17.75       -35.75       5 E017S36BB      NA       NA
#> 6        18.25       -35.75       6 E018S36BB      NA       NA
#>                         geometry             centroid
#> 1 POLYGON ((15.5 -36, 16 -36,... POINT (15.75 -35.75)
#> 2 POLYGON ((16 -36, 16.5 -36,... POINT (16.25 -35.75)
#> 3 POLYGON ((16.5 -36, 17 -36,... POINT (16.75 -35.75)
#> 4 POLYGON ((17 -36, 17.5 -36,... POINT (17.25 -35.75)
#> 5 POLYGON ((17.5 -36, 18 -36,... POINT (17.75 -35.75)
#> 6 POLYGON ((18 -36, 18.5 -36,... POINT (18.25 -35.75)
dim(grid_spp); head(grid_spp[, 1:8])
#> [1]  415 2874
#> # A tibble: 6 × 8
#>   grid_id centroid_lon centroid_lat mapsheet  obs_sum spp_rich
#>   <chr>          <dbl>        <dbl> <chr>       <dbl>    <dbl>
#> 1 1026            28.8        -22.3 E028S23BB       3        2
#> 2 1027            29.2        -22.3 E029S23BB      41       31
#> 3 1028            29.7        -22.3 E029S23BB      10       10
#> 4 1029            30.3        -22.3 E030S23BB       7        7
#> 5 1030            30.8        -22.3 E030S23BB       6        6
#> 6 1031            31.3        -22.3 E031S23BB     107       76
#> # ℹ 2 more variables: `Mylothris agathina subsp. agathina` <dbl>,
#> #   `Pieris brassicae` <dbl>
dim(grid_spp_pa); head(grid_spp_pa[, 1:8])
#> [1]  415 2874
#> # A tibble: 6 × 8
#>   grid_id centroid_lon centroid_lat mapsheet  obs_sum spp_rich
#>   <chr>          <dbl>        <dbl> <chr>       <dbl>    <dbl>
#> 1 1026            28.8        -22.3 E028S23BB       3        2
#> 2 1027            29.2        -22.3 E029S23BB      41       31
#> 3 1028            29.7        -22.3 E029S23BB      10       10
#> 4 1029            30.3        -22.3 E030S23BB       7        7
#> 5 1030            30.8        -22.3 E030S23BB       6        6
#> 6 1031            31.3        -22.3 E031S23BB     107       76
#> # ℹ 2 more variables: `Mylothris agathina subsp. agathina` <dbl>,
#> #   `Pieris brassicae` <dbl>
```

`grid_spp` now serves as the **site‑level backbone** for modelling
(e.g. spatial GLMs) or visualisation (e.g. dot plots), whereas
`grid_spp_pa` slots directly into Jaccard- or Sørensen-based
beta-diversity workflows. `site_spp` retains the raw observation detail
for drill‑down analyses.

#### 3. Visualise observation density across South Africa

With the grid summaries in hand we can now **map the spatial
distribution of observation effort**. The recipe below layers three
geometric objects in a single `ggplot2` call:

1.  **Grid polygons (`grid_sf`)**: Outlined in semi‑transparent grey to
    give a subtle sense of the analytical lattice without overwhelming
    the figure.
2.  **Centroid points (`grid_spp`)**: Plotted using longitude/latitude
    coordinates and symbol attributes that encode sampling intensity.
    For example, below **size & colour** are mapped to `sqrt(obs_sum)`.
    We use [`sqrt()`](https://rdrr.io/r/base/MathFun.html) because a
    square‑root transform is often preferable when counts span large
    orders of magnitude as it compresses large values while still
    highlighting structure among sparsely sampled cells.
3.  **National border (`rsa`)**: Emphasised in solid black to anchor the
    map in a familiar outline.

A perceptually uniform `Viridis` palette (`option = "turbo"`) supports
colour‑blind accessibility, while
[`theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
removes visual clutter so the data can speak for themselves.

``` r
ggplot2::ggplot() +
  # 1. grid polygons as subtle backdrop 
  geom_sf(data = grid_sf, fill = NA, colour = "darkgrey", linewidth = 0.2, alpha = 0.5) +
  
  # 2. centroids sized/coloured by sampling effort 
  geom_point(
    data = grid_spp,
    aes(x = centroid_lon, y = centroid_lat,
        size  = sqrt(obs_sum),
        colour = sqrt(obs_sum)),
    alpha = 0.8
  ) +
  
  # Divergent colour scale 
  scale_colour_viridis_c(option = "turbo", name = "√ Observations") +
  scale_size_continuous(name = "√ Observations", guide = "none") +
  
  # 3. national outline 
  geom_sf(data = rsa, fill = NA, colour = "black", linewidth = 0.4) +
  
  theme_minimal() +
  labs(
    title = "Observation density across South Africa (0.5° grid)",
    x = "Longitude", y = "Latitude"
  )
```

![Observation density across South Africa (0.5°
grid)](figures/2-map-aoi-1.png)

------------------------------------------------------------------------

#### 4. Visualise sampling effort and richness

[`generate_grid()`](https://b-cubed-eu.github.io/dissmapr/reference/generate_grid.md)
also returns a three-layer `SpatRaster` (`grid_r`) whose second and
third bands store cell-level metrics:

- `obs_sum`: Total observations aggregated across the chosen species
  columns (units = observation count)
- `spp_rich`: Number of species (non-zero columns) recorded in the cell
  (units = unique species count)

The chunk below extracts those two layers, applies a square-root stretch
(to dampen the influence of very large counts), and renders them
side-by-side with a perceptually uniform turbo palette.

``` r
# 1. Extract & transform layers (use terra method explicitly)
effRich_r = sqrt(grid_r[[c("obs_sum", "spp_rich")]])

# 2. Save and reset graphics state safely
old_par = par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

layout(matrix(1:2, nrow = 1))

titles = c(
  "Sampling effort (√obs count)",
  "Species richness (√unique count)"
)

for (i in 1:2) {
  terra::plot(
    effRich_r[[i]],
    col    = viridisLite::turbo(100),
    breaks = 100,
    colNA  = NA_character_,
    axes   = FALSE,
    main   = titles[i],
    cex.main = 0.8
  )
  terra::plot(terra::vect(rsa), add = TRUE, border = "black", lwd = 0.4)
}
```

![Sampling effort (√obs count) and Species richness (√unique
count)](figures/2-eff-rich-1.png)

These maps quickly reveal where sampling effort is concentrated and how
species richness varies across the landscape—useful diagnostics before
any downstream modelling.

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_South Africa.utf8  LC_CTYPE=English_South Africa.utf8   
#> [3] LC_MONETARY=English_South Africa.utf8 LC_NUMERIC=C                         
#> [5] LC_TIME=English_South Africa.utf8    
#> 
#> time zone: Africa/Johannesburg
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.3.0            pbapply_1.7-4        geodata_0.6-9       
#>   [4] pROC_1.19.0.1        s2_1.1.11            permute_0.9-10      
#>   [7] rlang_1.2.0          magrittr_2.0.5       otel_0.2.0          
#>  [10] e1071_1.7-17         compiler_4.5.2       mgcv_1.9-4          
#>  [13] systemfonts_1.3.2    vctrs_0.7.3          maps_3.4.3          
#>  [16] reshape2_1.4.5       stringr_1.6.0        wk_0.9.5            
#>  [19] pkgconfig_2.0.3      fastmap_1.2.0        labeling_0.4.3      
#>  [22] utf8_1.2.6           rmarkdown_2.31       prodlim_2026.03.11  
#>  [25] dissmapr_0.2.0       ragg_1.5.2           purrr_1.2.2         
#>  [28] xfun_0.59            cachem_1.1.0         jsonlite_2.0.0      
#>  [31] recipes_1.3.3        terra_1.9-34         parallel_4.5.2      
#>  [34] cluster_2.1.8.2      R6_2.6.1             bslib_0.11.0        
#>  [37] stringi_1.8.7        RColorBrewer_1.1-3   parallelly_1.47.0   
#>  [40] rpart_4.1.27         lubridate_1.9.5      jquerylib_0.1.4     
#>  [43] estimability_1.5.1   Rcpp_1.1.1-1.1       iterators_1.0.14    
#>  [46] knitr_1.51           future.apply_1.20.2  fields_17.3         
#>  [49] zoo_1.8-15           Matrix_1.7-5         splines_4.5.2       
#>  [52] nnet_7.3-20          timechange_0.4.0     tidyselect_1.2.1    
#>  [55] rstudioapi_0.19.0    yaml_2.3.12          vegan_2.7-5         
#>  [58] timeDate_4052.112    codetools_0.2-20     listenv_1.0.0       
#>  [61] lattice_0.22-9       tibble_3.3.1         plyr_1.8.9          
#>  [64] withr_3.0.3          S7_0.2.2             geosphere_1.6-8     
#>  [67] evaluate_1.0.5       sf_1.1-1             future_1.70.0       
#>  [70] desc_1.4.3           survival_3.8-6       units_1.0-1         
#>  [73] proxy_0.4-29         mclust_6.1.2         pillar_1.11.1       
#>  [76] KernSmooth_2.23-26   corrplot_0.95        renv_1.1.4          
#>  [79] foreach_1.5.2        stats4_4.5.2         generics_0.1.4      
#>  [82] zetadiv_1.3.0        scales_1.4.0         globals_0.19.1      
#>  [85] xtable_1.8-8         class_7.3-23         glue_1.8.1          
#>  [88] clValid_0.7          emmeans_2.0.3        tools_4.5.2         
#>  [91] data.table_1.18.4    ModelMetrics_1.2.2.2 gower_1.0.2         
#>  [94] fs_2.1.0             mvtnorm_1.4-1        dotCall64_1.2       
#>  [97] grid_4.5.2           tidyr_1.3.2          ipred_0.9-15        
#> [100] patchwork_1.3.2      nlme_3.1-169         cli_3.6.6           
#> [103] rappdirs_0.3.4       textshaping_1.0.5    NbClust_3.0.1       
#> [106] spam_2.11-4          viridisLite_0.4.3    lava_1.9.1          
#> [109] scam_1.2-22          dplyr_1.2.1          gtable_0.3.6        
#> [112] sass_0.4.10          digest_0.6.39        classInt_0.4-11     
#> [115] caret_7.0-1          ggrepel_0.9.8        htmlwidgets_1.6.4   
#> [118] farver_2.1.2         entropy_1.3.2        htmltools_0.5.9     
#> [121] pkgdown_2.2.0        lifecycle_1.0.5      factoextra_2.0.0    
#> [124] httr_1.4.8           hardhat_1.4.3        MASS_7.3-65
```
