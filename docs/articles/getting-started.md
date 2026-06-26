# getting-started

------------------------------------------------------------------------

## `dissmapr`

### A Novel Framework for Automated Compositional Dissimilarity and Biodiversity Turnover Analysis

------------------------------------------------------------------------

### Introduction

`dissmapr` is an R package for analysing compositional dissimilarity and
biodiversity turnover across spatial gradients. It provides scalable,
modular workflows that integrate species occurrence, environmental data,
and multi-site compositional turnover metrics to quantify and predict
biodiversity patterns. A core feature is the use of zeta diversity,
which extends beyond pairwise comparisons to capture shared species
across multiple sites - offering deeper insight into community assembly,
turnover, and connectivity, for both rare and common species. By
incorporating different regression methods within the framework of
Multi-Site Generalised Dissimilarity Modelling (MS-GDM), `dissmapr`
enables robust mapping, bioregional classification, and scenario-based
forecasting. Designed for flexibility and reproducibility, it supports
biodiversity monitoring and conservation planning at landscape to
regional scales.

------------------------------------------------------------------------

### Step-by-Step Workflow

`dissmapr` implements a structured, reproducible workflow for analysing
biodiversity patterns and delineating bioregions. Each function aligns
with a specific step in the workflow, guiding users from data
acquisition to predictive mapping. The workflow begins with sourcing
species occurrence and georeferenced environmental data, followed by
data formatting and calculation of compositional turnover using zeta
diversity metrics (via the `zetadiv` package). Multi-Site Generalized
Dissimilarity Modelling (MS-GDM) is then applied to model and predict
dissimilarity across landscapes. As part of the Dissimilarity Cube in
the [Biodiversity Building Blocks for Policy](https://b-cubed.eu/)
project, the outputs and predictions from MS-GDM are classified into
spatial clusters of species composition, representing distinct
bioregions. The workflow also supports the integration of historical and
future climate data to assess shifts in biodiversity and detect
emerging, shifting, or dissolving bioregions under global change. This
step-by-step structure, mirrored in the accompanying tutorial sections,
promotes accessibility, transparency, and ecological insight at multiple
spatial and temporal scales.

------------------------------------------------------------------------

#### Install and load `dissmapr`

Install and load the `dissmapr` package from GitHub, ensuring all
functions are available for use in the workflow.

``` r
# install remotes if needed
# install.packages("remotes")
remotes::install_github("macSands/dissmapr")
```

``` r
# Ensure the package is loaded when knitting
library(dissmapr)

# Make sure all the functions are loaded
devtools::load_all()
```

------------------------------------------------------------------------

#### Load other R libraries

Load core libraries for spatial processing, biodiversity modelling, and
visualization required across the `dissmapr` analysis pipeline.

``` r
# Load necessary libraries
library(httr)       # HTTP client  
library(geodata)    # Download geographic data  
library(data.table) # Fast large-table operations  
library(dplyr)      # Data manipulation verbs  
library(tidyr)      # Tidy data reshaping  
library(zoo)        # Time series utilities  
library(sf)         # Vector spatial data  
library(terra)      # Raster spatial operations  
library(tidyterra)  # supplies geom_spatraster()
library(zetadiv)    # Multi-site dissimilarity modelling
library(ggplot2)    # Grammar of graphics  
library(viridis)    # Perceptual color scales  
```

------------------------------------------------------------------------

#### Get species occurrence records using `get_occurrence_data()`

To contextualise the following steps of the workflow, we use the
butterfly data in South Africa, accessed from GBIF, as the case for
demonstration. Ultimately, the choice for the AoI and taxa is
user-specific. This section focuses on automating the retrieval and
pre-processing of biodiversity occurrence data from various sources,
including:

1.  local `databases`.csv\` files,
2.  URLs or `.zip` files from the Global Biodiversity Information
    Facility (GBIF), and
3.  species occurrence cubes from B3 (specification) \[*work in
    progress*\].

The function assembles data on species distributions across specified
taxonomic groups and regions, producing presence-absence or abundance
matrices that quantify species co-occurrence within locations.

``` r
bfly_data = get_occurrence_data(
  data        = system.file("extdata", "gbif_butterflies.csv", package = "dissmapr"),
  source_type = 'local_csv',
  sep         = '\t'
)

# Check results
dim(bfly_data)
#> [1] 81825    52
str(bfly_data)
#> 'data.frame':    81825 obs. of  52 variables:
#>  $ gbifID                          : num  9.23e+08 9.23e+08 9.23e+08 9.22e+08 9.22e+08 ...
#>  $ datasetKey                      : chr  "6ac3f774-d9fb-4796-b3e9-92bf6c81c084" "6ac3f774-d9fb-4796-b3e9-92bf6c81c084" "6ac3f774-d9fb-4796-b3e9-92bf6c81c084" "6ac3f774-d9fb-4796-b3e9-92bf6c81c084" ...
#>  $ occurrenceID                    : chr  "" "" "" "" ...
#>  $ kingdom                         : chr  "Animalia" "Animalia" "Animalia" "Animalia" ...
#>  $ phylum                          : chr  "Arthropoda" "Arthropoda" "Arthropoda" "Arthropoda" ...
#>  $ class                           : chr  "Insecta" "Insecta" "Insecta" "Insecta" ...
#>  $ order                           : chr  "Lepidoptera" "Lepidoptera" "Lepidoptera" "Lepidoptera" ...
#>  $ family                          : chr  "Pieridae" "Pieridae" "Papilionidae" "Pieridae" ...
#>  $ genus                           : chr  "Pieris" "Pieris" "Papilio" "Mylothris" ...
#>  $ sp_name                         : chr  "Pieris brassicae" "Pieris brassicae" "Papilio demodocus" "Mylothris agathina" ...
#>  $ infraspecificEpithet            : chr  "" "" "" "agathina" ...
#>  $ taxonRank                       : chr  "SPECIES" "SPECIES" "SPECIES" "SUBSPECIES" ...
#>  $ scientificName                  : chr  "Pieris brassicae (Linnaeus, 1758)" "Pieris brassicae (Linnaeus, 1758)" "Papilio demodocus Esper, 1798" "Mylothris agathina agathina" ...
#>  $ verbatimScientificName          : chr  "Pieris brassicae" "Pieris brassicae" "Papilio demodocus subsp. demodocus" "Mylothris agathina subsp. agathina" ...
#>  $ verbatimScientificNameAuthorship: logi  NA NA NA NA NA NA ...
#>  $ countryCode                     : chr  "ZA" "ZA" "ZA" "ZA" ...
#>  $ locality                        : chr  "Hermanus" "Polkadraai Road" "Signal Hill" "Hermanus" ...
#>  $ stateProvince                   : chr  "" "" "" "" ...
#>  $ occurrenceStatus                : chr  "PRESENT" "PRESENT" "PRESENT" "PRESENT" ...
#>  $ individualCount                 : int  NA NA NA NA NA NA NA NA NA NA ...
#>  $ publishingOrgKey                : chr  "bb646dff-a905-4403-a49b-6d378c2cf0d9" "bb646dff-a905-4403-a49b-6d378c2cf0d9" "bb646dff-a905-4403-a49b-6d378c2cf0d9" "bb646dff-a905-4403-a49b-6d378c2cf0d9" ...
#>  $ y                               : num  -34.4 -34 -33.9 -34.4 -34.4 ...
#>  $ x                               : num  19.2 18.8 18.4 19.2 18.5 ...
#>  $ coordinateUncertaintyInMeters   : num  250 250 250 250 250 250 250 250 250 250 ...
#>  $ coordinatePrecision             : logi  NA NA NA NA NA NA ...
#>  $ elevation                       : logi  NA NA NA NA NA NA ...
#>  $ elevationAccuracy               : logi  NA NA NA NA NA NA ...
#>  $ depth                           : logi  NA NA NA NA NA NA ...
#>  $ depthAccuracy                   : logi  NA NA NA NA NA NA ...
#>  $ eventDate                       : chr  "2012-10-13T00:00" "2012-11-01T00:00" "2012-10-31T00:00" "2012-10-13T00:00" ...
#>  $ day                             : int  13 1 31 13 30 23 20 22 21 22 ...
#>  $ month                           : int  10 11 10 10 10 10 10 10 10 10 ...
#>  $ year                            : int  2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 ...
#>  $ taxonKey                        : int  1920506 1920506 1938125 11374894 6113873 6118817 6118817 6118871 1919924 1770649 ...
#>  $ speciesKey                      : int  1920506 1920506 1938125 5137998 6113873 1954382 1954382 1964076 1919924 1770649 ...
#>  $ basisOfRecord                   : chr  "HUMAN_OBSERVATION" "HUMAN_OBSERVATION" "HUMAN_OBSERVATION" "HUMAN_OBSERVATION" ...
#>  $ institutionCode                 : chr  "naturgucker" "naturgucker" "naturgucker" "naturgucker" ...
#>  $ collectionCode                  : chr  "naturgucker" "naturgucker" "naturgucker" "naturgucker" ...
#>  $ catalogNumber                   : chr  "-190723441" "2051102952" "1565945073" "-1879976251" ...
#>  $ recordNumber                    : chr  "" "" "" "" ...
#>  $ identifiedBy                    : chr  "" "" "" "" ...
#>  $ dateIdentified                  : chr  "" "" "" "" ...
#>  $ license                         : chr  "CC_BY_4_0" "CC_BY_4_0" "CC_BY_4_0" "CC_BY_4_0" ...
#>  $ rightsHolder                    : chr  "" "" "" "" ...
#>  $ recordedBy                      : chr  "591374253" "591374253" "591374253" "591374253" ...
#>  $ typeStatus                      : logi  NA NA NA NA NA NA ...
#>  $ establishmentMeans              : logi  NA NA NA NA NA NA ...
#>  $ lastInterpreted                 : chr  "2024-03-15T23:20:20.346Z" "2024-03-15T23:22:38.206Z" "2024-03-15T23:26:59.721Z" "2024-03-15T23:23:55.839Z" ...
#>  $ mediaType                       : chr  "" "" "" "StillImage" ...
#>  $ issue                           : chr  "COORDINATE_ROUNDED;GEODETIC_DATUM_ASSUMED_WGS84;CONTINENT_DERIVED_FROM_COORDINATES;MULTIMEDIA_URI_INVALID" "COORDINATE_ROUNDED;GEODETIC_DATUM_ASSUMED_WGS84;CONTINENT_DERIVED_FROM_COORDINATES;MULTIMEDIA_URI_INVALID" "COORDINATE_ROUNDED;GEODETIC_DATUM_ASSUMED_WGS84;CONTINENT_DERIVED_FROM_COORDINATES;TAXON_MATCH_HIGHERRANK;MULTI"| __truncated__ "COORDINATE_ROUNDED;GEODETIC_DATUM_ASSUMED_WGS84;CONTINENT_DERIVED_FROM_COORDINATES" ...
#>  $ site_id                         : int  1 2 3 1 4 5 5 5 5 5 ...
#>  $ pa                              : num  1 1 1 1 1 1 1 1 1 1 ...
head(bfly_data[,1:6])
#>      gbifID                           datasetKey occurrenceID  kingdom
#> 1 923051749 6ac3f774-d9fb-4796-b3e9-92bf6c81c084              Animalia
#> 2 922985630 6ac3f774-d9fb-4796-b3e9-92bf6c81c084              Animalia
#> 3 922619348 6ac3f774-d9fb-4796-b3e9-92bf6c81c084              Animalia
#> 4 922426210 6ac3f774-d9fb-4796-b3e9-92bf6c81c084              Animalia
#> 5 921650584 6ac3f774-d9fb-4796-b3e9-92bf6c81c084              Animalia
#> 6 921485695 6ac3f774-d9fb-4796-b3e9-92bf6c81c084              Animalia
#>       phylum   class
#> 1 Arthropoda Insecta
#> 2 Arthropoda Insecta
#> 3 Arthropoda Insecta
#> 4 Arthropoda Insecta
#> 5 Arthropoda Insecta
#> 6 Arthropoda Insecta
```

------------------------------------------------------------------------

#### Format data using `format_df()`

Use `format_df` to *standardise and reshape* raw biodiversity tables
into the *long* or *wide* format required by later `dissmapr` steps.
Importantly, this function does not alter the spatial resolution of the
original observations - it simply tidies the data by automatically
identifying key columns (e.g., coordinates, species, and values),
assigning missing observation/site IDs, and reformatting the data for
analysis. Outputs include a cleaned `site_obs` dataset and `site_spp`
matrix for further processing:

• **site_obs**: Simplified table with unique `site_id`, `x`, `y`,
`species` and `value` records (long format). • **site_spp**:
Site-by-species matrix for biodiversity assessments (wide format).

##### Format data into long (site_obs) and wide (site_spp) formats

``` r
bfly_result = format_df(
  data        = bfly_data, # A `data.frame` of biodiversity records
  species_col = 'verbatimScientificName', # Name of species column (required for `"long"`)
  value_col   = 'pa', # Name of value column (e.g. presence/abundance; for `"long"`)
  extra_cols  = NULL, # Character vector of other columns to keep
  format      = 'long' # Either`"long"` or `"wide"`. If `NULL`, inferred from `species_col` & `value_col`
)

# Check `bfly_result` structure
str(bfly_result, max.level = 1)
#> List of 2
#>  $ site_obs:'data.frame':    79953 obs. of  5 variables:
#>  $ site_spp: tibble [56,090 × 2,871] (S3: tbl_df/tbl/data.frame)

# Optional: Create new objects from list items
site_obs = bfly_result$site_obs
site_spp = bfly_result$site_spp

# Check results
dim(site_obs)
#> [1] 79953     5
head(site_obs)
#>   site_id        x         y                            species value
#> 1       1 19.24410 -34.42086                   Pieris brassicae     1
#> 2       2 18.75564 -33.96044                   Pieris brassicae     1
#> 3       3 18.40321 -33.91651 Papilio demodocus subsp. demodocus     1
#> 4       1 19.24410 -34.42086 Mylothris agathina subsp. agathina     1
#> 5       4 18.47488 -34.35024                  Eutricha capensis     1
#> 6       5 25.65097 -33.58570            Drepanogynis bifasciata     1

dim(site_spp)
#> [1] 56090  2871
head(site_spp[,1:6])
#> # A tibble: 6 × 6
#>   site_id     x     y `Mylothris agathina subsp. agathina` `Pieris brassicae`
#>     <int> <dbl> <dbl>                                <dbl>              <dbl>
#> 1       1  19.2 -34.4                                    1                  1
#> 2       2  18.8 -34.0                                    0                  1
#> 3       3  18.4 -33.9                                    0                  0
#> 4       4  18.5 -34.4                                    0                  0
#> 5       5  25.7 -33.6                                    0                  0
#> 6       6  22.2 -33.6                                    0                  0
#> # ℹ 1 more variable: `Tarucus thespis` <dbl>

#### Get parameters from processed data to use later
# Number of species
(n_sp = dim(site_spp)[2] - 3)
#> [1] 2868

# Species names
sp_cols = names(site_spp)[-c(1:3)]
sp_cols[1:10]
#>  [1] "Mylothris agathina subsp. agathina" "Pieris brassicae"                  
#>  [3] "Tarucus thespis"                    "Acraea horta"                      
#>  [5] "Danaus chrysippus"                  "Papilio demodocus subsp. demodocus"
#>  [7] "Eutricha capensis"                  "Mesocelis monticola"               
#>  [9] "Vanessa cardui"                     "Cuneisigna obstans"
```

------------------------------------------------------------------------

#### User-defined area of interest and grid resolution

Load the spatial boundary data for South Africa to serve as the
geographic reference for all subsequent biodiversity analyses and
visualizations.

``` r
# Read RSA shape file
rsa = sf::st_read(system.file("extdata", "rsa.shp", package = "dissmapr"))
#> Reading layer `rsa' from data source 
#>   `D:\Methods\R\myR_Packages\myCompletePks\dissmapr\inst\extdata\rsa.shp' 
#>   using driver `ESRI Shapefile'
#> Simple feature collection with 1 feature and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 16.45802 ymin: -34.83514 xmax: 32.89125 ymax: -22.12661
#> Geodetic CRS:  WGS 84

# Define your resolution and create mask to use later
res = 0.5 # 0.5 degrees is roughly 55km

# Convert to a terra vector
rsa_vect = vect(rsa)

# Create an empty raster over RSA at your desired resolution
grid = rast(rsa_vect, resolution = res, crs = crs(rsa_vect))
values(grid) = 1   # fill with dummy values

# Mask everything outside the RSA boundary
grid_masked = mask(grid, rsa_vect)
```

------------------------------------------------------------------------

#### Summarise records by grid centroid using `generate_grid()`

Use `generate_grid` to overlay a user-defined lattice on the study
region—whose bounds are inferred from the occurrence data—without
modifying the underlying observations themselves. The function
constructs a raster‐based grid of the chosen cell size, assigns each
point a unique `grid_id`, and compiles a summary of user-specified
attributes for every cell. It returns three coordinated outputs:

• **grid**: `SpatRaster` with grid index • **grid_sf**: `sf` and
`data.frame` i.e. lattice polygon features for mapping or spatial joins
• **block_sp**: `data.frame` that records per-cell totals, centroids,
and other statistics.

By aggregating raw records into consistent spatial units,
`generate_grid` provides the structured foundation needed for subsequent
landscape-scale biodiversity analyses.

------------------------------------------------------------------------

##### Assign records to a grid at a set resolution

Aggregate species records into grid cells of user-specified size
(e.g. 0.5°) to enable spatially standardized biodiversity analyses.

``` r
grid_list = generate_grid(
  data       = site_spp,
  x_col      = "x",
  y_col      = "y",
  grid_size  = 0.5,
  sum_col_range = 4:ncol(site_spp),
  crs_epsg   = 4326
)

# Check `grid_list` structure
str(grid_list, max.level = 1)
#> List of 3
#>  $ grid    :S4 class 'SpatRaster' [package "terra"]
#>  $ grid_sf :Classes 'sf' and 'data.frame':   1110 obs. of  8 variables:
#>   ..- attr(*, "sf_column")= chr "geometry"
#>   ..- attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA NA NA NA NA NA
#>   .. ..- attr(*, "names")= chr [1:7] "centroid_lon" "centroid_lat" "grid_id" "mapsheet" ...
#>  $ block_sp:'data.frame':    415 obs. of  2874 variables:

# Optional: Create new objects from list items
aoi_grid = grid_list$grid_sf
grid_spp = grid_list$block_sp

# Check results
dim(aoi_grid)
#> [1] 1110    8
head(aoi_grid)
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

dim(grid_spp)
#> [1]  415 2874
head(grid_spp[,1:6])
#>   grid_id centroid_lon centroid_lat  mapsheet obs_sum spp_rich
#> 1    1026        28.75    -22.25004 E028S23BB       3        2
#> 2    1027        29.25    -22.25004 E029S23BB      41       31
#> 3    1028        29.75    -22.25004 E029S23BB      10       10
#> 4    1029        30.25    -22.25004 E030S23BB       7        7
#> 5    1030        30.75    -22.25004 E030S23BB       6        6
#> 6    1031        31.25    -22.25004 E031S23BB     107       76
```

##### Generate a data frame ‘xy’ of site centroids

Extract longitude–latitude coordinates and summary metrics for each
occupied grid cell.

``` r
# Grid centroids with 'gird_id', 'centroid_lon', 'centroid_lat', 'obs_sum' and `spp_rich`
grid_xy = grid_spp[,c(1:3,5:6)]

# Create species observations data.frame
spp_obs = site_obs

# Check results
dim(grid_xy)
#> [1] 415   5
head(grid_xy)
#>   grid_id centroid_lon centroid_lat obs_sum spp_rich
#> 1    1026        28.75    -22.25004       3        2
#> 2    1027        29.25    -22.25004      41       31
#> 3    1028        29.75    -22.25004      10       10
#> 4    1029        30.25    -22.25004       7        7
#> 5    1030        30.75    -22.25004       6        6
#> 6    1031        31.25    -22.25004     107       76

dim(spp_obs)
#> [1] 79953     5
head(spp_obs)
#>   site_id        x         y                            species value
#> 1       1 19.24410 -34.42086                   Pieris brassicae     1
#> 2       2 18.75564 -33.96044                   Pieris brassicae     1
#> 3       3 18.40321 -33.91651 Papilio demodocus subsp. demodocus     1
#> 4       1 19.24410 -34.42086 Mylothris agathina subsp. agathina     1
#> 5       4 18.47488 -34.35024                  Eutricha capensis     1
#> 6       5 25.65097 -33.58570            Drepanogynis bifasciata     1
```

##### Generate a map of RSA with occupied grid cells as centroid points

Visualize observation density across South Africa using centroid-based
mapping of gridded biodiversity data.

``` r
ggplot() +
  geom_sf(data = aoi_grid, fill = NA, color = "darkgrey", alpha = 0.5) +
  geom_point(data = grid_spp,
             aes(x = centroid_lon, y = centroid_lat,
                 size = sqrt(obs_sum),
                 color = sqrt(obs_sum))) +
  scale_color_viridis_c(option = "turbo") +
  geom_sf(data = rsa, fill = NA, color = "black") +
  theme_minimal() +
  labs(title = "0.5° Grid with Observation Counts",
       x = "Longitude", y = "Latitude")
```

![](reference/figures/full-map-aoi-1.png)

------------------------------------------------------------------------
