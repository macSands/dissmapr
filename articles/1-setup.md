# Getting started

## `dissmapr`

### Framework for Compositional Dissimilarity and Biodiversity Turnover Analysis

#### Introduction

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

#### 1. Install and load `dissmapr`

Install and load the `dissmapr` package from GitHub, ensuring all
functions are available for use in the workflow.

``` r

# install remotes if needed
# install.packages("remotes")
# remotes::install_github("macSands/dissmapr")
```

``` r

# Ensure the package is loaded when knitting
library(dissmapr)

# Make sure all the functions are loaded
# devtools::load_all()
```

#### 2. Load other R libraries

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
library(patchwork)  # Sequentially build up plots on one page
library(mclust)     # Clustering, Classification, and Density Estimation
```

#### 3. Get species occurrence records using `get_occurrence_data()`

To contextualise the following steps of the workflow, we use South
African butterfly data accessed from GBIF ([DOI:
10.15468/dl.jh6maj](https://www.gbif.org/occurrence/download/0006880-241024112534372)),
as a demonstration case. Ultimately, the choice for the Area of Interest
(AoI) and taxa is user-specific. This section demonstrates how to
automate the retrieval and pre-processing of biodiversity occurrence
data from a GBIF query (stored locally as a `.csv` file), however the
same workflow can ingest other sources as well (see the
[`get_occurrence_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_occurrence_data.md)
documentation for details). Data inputs currently supported include:

- **Local** databases or `.csv` files
- **URLs** or `.zip` files from the Global Biodiversity Information
  Facility (GBIF)
- Future inclusion of **GBIF species occurrence cubes**. Read the
  [species occurrence cubes in
  GBIF](https://www.gbif.org/occurrence-cubes) documentation for full
  details on creating, customizing and submitting queries for occurrence
  cubes. Read the [b-cubed](https://b-cubed.eu/) documentation on
  [specification for species occurrence cubes and their
  production](https://docs.b-cubed.eu/guides/occurrence-cube/).

[`get_occurrence_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_occurrence_data.md)
then organises the records by the chosen taxonomic scope and region,
returning presence–absence and/or abundance matrices that summarise
species co-occurrence records with latitude and longitude coordinates.

``` r

load(system.file("extdata", "gbif_butterflies_csv.RData", package = "dissmapr"), envir = knitr::knit_global())

bfly_data = get_occurrence_data(
  data        = gbif_butterflies_csv,
  source_type = 'data_frame'
)

# bfly_data = get_occurrence_data(
#   data        = system.file("extdata", "gbif_butterflies.csv", package = "dissmapr"),
#   source_type = 'local_csv',
#   sep         = '\t'
# )

# Check results but only a subset of columns to fit in console
dim(bfly_data)
#> [1] 81825    52
str(bfly_data[,c(51,52,22,23,1,14,16,17,30)]) 
#> 'data.frame':    81825 obs. of  9 variables:
#>  $ site_id               : int  1 2 3 1 4 5 5 5 5 5 ...
#>  $ pa                    : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ y                     : num  -34.4 -34 -33.9 -34.4 -34.4 ...
#>  $ x                     : num  19.2 18.8 18.4 19.2 18.5 ...
#>  $ gbifID                : num  9.23e+08 9.23e+08 9.23e+08 9.22e+08 9.22e+08 ...
#>  $ verbatimScientificName: chr  "Pieris brassicae" "Pieris brassicae" "Papilio demodocus subsp. demodocus" "Mylothris agathina subsp. agathina" ...
#>  $ countryCode           : chr  "ZA" "ZA" "ZA" "ZA" ...
#>  $ locality              : chr  "Hermanus" "Polkadraai Road" "Signal Hill" "Hermanus" ...
#>  $ eventDate             : chr  "2012-10-13T00:00" "2012-11-01T00:00" "2012-10-31T00:00" "2012-10-13T00:00" ...
head(bfly_data[,c(51,52,22,23,1,14,16,17,30)])
#>   site_id pa         y        x    gbifID             verbatimScientificName
#> 1       1  1 -34.42086 19.24410 923051749                   Pieris brassicae
#> 2       2  1 -33.96044 18.75564 922985630                   Pieris brassicae
#> 3       3  1 -33.91651 18.40321 922619348 Papilio demodocus subsp. demodocus
#> 4       1  1 -34.42086 19.24410 922426210 Mylothris agathina subsp. agathina
#> 5       4  1 -34.35024 18.47488 921650584                  Eutricha capensis
#> 6       5  1 -33.58570 25.65097 921485695            Drepanogynis bifasciata
#>   countryCode                                          locality
#> 1          ZA                                          Hermanus
#> 2          ZA                                   Polkadraai Road
#> 3          ZA                                       Signal Hill
#> 4          ZA                                          Hermanus
#> 5          ZA Cape of Good Hope / Cape Point Area, South Africa
#> 6          ZA                             Kudu Ridge Game Lodge
#>          eventDate
#> 1 2012-10-13T00:00
#> 2 2012-11-01T00:00
#> 3 2012-10-31T00:00
#> 4 2012-10-13T00:00
#> 5 2012-10-30T00:00
#> 6 2012-10-23T00:00
```

#### 4. Format data using `format_df()`

Use
[`format_df()`](https://b-cubed-eu.github.io/dissmapr/reference/format_df.md)
to *standardise and reshape* raw biodiversity tables into the *long* or
*wide* format required by later `dissmapr` steps. Importantly, this
function does not alter the spatial resolution of the original
observations - it simply tidies the data by automatically identifying
key columns (e.g., coordinates, species, and values), assigning unique
site IDs (`site_id`), renaming or removing columns, and reformatting the
data for analysis. Outputs include a cleaned `site_obs` dataset and
`site_spp` matrix for further processing:

- **site_obs**: Simplified table with unique `site_id`, `x`, `y`,
  `species` and `value` records (long format).
- **site_spp**: Site-by-species matrix for biodiversity assessments
  (wide format).

**Format data into long (`site_obs`) and wide (`site_spp`) formats**

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
#>  $ site_spp:'data.frame':    56090 obs. of  2871 variables:

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
#>   site_id        x         y Mylothris agathina subsp. agathina
#> 1       1 19.24410 -34.42086                                  1
#> 2       2 18.75564 -33.96044                                  0
#> 3       3 18.40321 -33.91651                                  0
#> 4       4 18.47488 -34.35024                                  0
#> 5       5 25.65097 -33.58570                                  0
#> 6       6 22.20197 -33.59240                                  0
#>   Pieris brassicae Tarucus thespis
#> 1                1               1
#> 2                1               0
#> 3                0               0
#> 4                0               0
#> 5                0               0
#> 6                0               0

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

``` r

# Each article in this series is self-contained: it loads the objects it needs
# from the single bundled snapshot `_data/dissmapr_vignettes.rds`, so the
# articles no longer have to be knitted in sequence.
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> other attached packages:
#>  [1] mclust_6.1.3      patchwork_1.3.2   viridis_0.6.5     viridisLite_0.4.3
#>  [5] ggplot2_4.0.3     zetadiv_1.3.0     scam_1.2-22       tidyterra_1.2.0  
#>  [9] sf_1.1-1          zoo_1.8-15        tidyr_1.3.2       dplyr_1.2.1      
#> [13] data.table_1.18.4 geodata_0.6-9     terra_1.9-34      httr_1.4.8       
#> [17] dissmapr_0.2.0   
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.3.0            pbapply_1.7-4        pROC_1.19.0.1       
#>   [4] gridExtra_2.3.1      permute_0.9-10       rlang_1.3.0         
#>   [7] magrittr_2.0.5       otel_0.2.0           e1071_1.7-17        
#>  [10] compiler_4.6.1       mgcv_1.9-4           systemfonts_1.3.2   
#>  [13] vctrs_0.7.3          maps_3.4.3           reshape2_1.4.5      
#>  [16] stringr_1.6.0        pkgconfig_2.0.3      fastmap_1.2.0       
#>  [19] rmarkdown_2.31       prodlim_2026.03.11   ragg_1.5.2          
#>  [22] purrr_1.2.2          xfun_0.60            cachem_1.1.0        
#>  [25] jsonlite_2.0.0       recipes_1.3.3        parallel_4.6.1      
#>  [28] cluster_2.1.8.2      R6_2.6.1             bslib_0.11.0        
#>  [31] stringi_1.8.7        RColorBrewer_1.1-3   parallelly_1.48.0   
#>  [34] rpart_4.1.27         estimability_2.0.0   lubridate_1.9.5     
#>  [37] jquerylib_0.1.4      Rcpp_1.1.2           iterators_1.0.14    
#>  [40] knitr_1.51           future.apply_1.20.2  fields_17.3         
#>  [43] Matrix_1.7-5         splines_4.6.1        nnet_7.3-20         
#>  [46] timechange_0.4.0     tidyselect_1.2.1     yaml_2.3.12         
#>  [49] vegan_2.7-5          timeDate_4052.112    codetools_0.2-20    
#>  [52] listenv_1.0.0        lattice_0.22-9       tibble_3.3.1        
#>  [55] plyr_1.8.9           withr_3.0.3          S7_0.2.2            
#>  [58] geosphere_1.6-8      evaluate_1.0.5       future_1.70.0       
#>  [61] desc_1.4.3           survival_3.8-6       units_1.0-1         
#>  [64] proxy_0.4-29         pillar_1.11.1        KernSmooth_2.23-26  
#>  [67] corrplot_0.95        renv_1.1.4           foreach_1.5.2       
#>  [70] stats4_4.6.1         generics_0.1.4       scales_1.4.0        
#>  [73] xtable_1.8-8         globals_0.19.1       class_7.3-23        
#>  [76] glue_1.8.1           clValid_0.7          emmeans_2.0.3       
#>  [79] tools_4.6.1          ModelMetrics_1.2.2.2 gower_1.0.2         
#>  [82] mvtnorm_1.4-2        fs_2.1.0             dotCall64_1.2       
#>  [85] grid_4.6.1           ipred_0.9-15         nlme_3.1-169        
#>  [88] cli_3.6.6            rappdirs_0.3.4       textshaping_1.0.5   
#>  [91] NbClust_3.0.1        spam_2.11-4          lava_1.9.2          
#>  [94] gtable_0.3.6         sass_0.4.10          digest_0.6.39       
#>  [97] classInt_0.4-11      caret_7.0-1          ggrepel_0.9.8       
#> [100] htmlwidgets_1.6.4    farver_2.1.2         entropy_1.3.2       
#> [103] htmltools_0.5.9      pkgdown_2.2.1        lifecycle_1.0.5     
#> [106] factoextra_2.1.0     hardhat_1.4.3        MASS_7.3-65
```
