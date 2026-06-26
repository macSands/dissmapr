# Package index

## Package

- [`dissmapr`](https://b-cubed-eu.github.io/dissmapr/reference/dissmapr-package.md)
  [`dissmapr-package`](https://b-cubed-eu.github.io/dissmapr/reference/dissmapr-package.md)
  : dissmapr: Workflow for Compositional Dissimilarity and Biodiversity
  Turnover Analysis

## Overview of dissmapr functions

The package consists of 10 core functions that address key steps in
biodiversity data analysis, from sourcing raw occurrence and
environmental data to mapping bioregions and predicting compositional
turnover. It also leverages existing tools, such as the `zetadiv`
package (specifically the
[`zetadiv::Zeta.msgdm`](https://rdrr.io/pkg/zetadiv/man/Zeta.msgdm.html)
function) for MS-GDM, alongside widely used R libraries like `terra`,
`sf`, and `ggplot2`.

## Core Functions

- [`get_occurrence_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_occurrence_data.md)
  : Import and harmonise biodiversity-occurrence data
- [`generate_grid()`](https://b-cubed-eu.github.io/dissmapr/reference/generate_grid.md)
  : Generate Spatial Grid and Gridded Summaries
- [`assign_mapsheet()`](https://b-cubed-eu.github.io/dissmapr/reference/assign_mapsheet.md)
  : Add Nearest Mapsheet Code and Center Coordinates
- [`get_enviro_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_enviro_data.md)
  : Retrieve, crop, resample, and link environmental rasters to sampling
  sites
- [`format_df()`](https://b-cubed-eu.github.io/dissmapr/reference/format_df.md)
  : Format Biodiversity Records to Long / Wide
- [`compute_orderwise()`](https://b-cubed-eu.github.io/dissmapr/reference/compute_orderwise.md)
  : Compute Order-wise Metrics
- [`rm_correlated()`](https://b-cubed-eu.github.io/dissmapr/reference/rm_correlated.md)
  : Remove Highly Correlated Predictors
- [`predict_dissim()`](https://b-cubed-eu.github.io/dissmapr/reference/predict_dissim.md)
  : Predict Pairwise Compositional Turnover (zeta-dissimilarity) with
  Richness
- [`map_bioreg()`](https://b-cubed-eu.github.io/dissmapr/reference/map_bioreg.md)
  : Raster-based Clustering and Interpolation of Bioregional Data
- [`map_bioregDiff()`](https://b-cubed-eu.github.io/dissmapr/reference/map_bioregDiff.md)
  : Map Bioregional Change Metrics Between Categorical Raster Layers

## ζ-MSGDM Workflow

Functions to automate a ζ-MSGDM workflow that fits, extracts and
visualizes I-spline models for any set of zeta orders in just three
function calls:

- [`run_ispline_models()`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md)
  : Run multiple Zeta.msgdm ispline models and return both models and
  combined ispline table
- [`plot_ispline_lines()`](https://b-cubed-eu.github.io/dissmapr/reference/plot_ispline_lines.md)
  : Plot ispline partial effects with quantile and start-point markers
- [`plot_ispline_boxplots()`](https://b-cubed-eu.github.io/dissmapr/reference/plot_ispline_boxplots.md)
  : Plot facetted boxplots for all ispline basis columns

## Available Metrics

Use `helper_indices` to choose a specific function for the func
parameter in `compute_orderwise`, enabling the calculation of one of the
following metrics:

- [`geodist_helper()`](https://b-cubed-eu.github.io/dissmapr/reference/geodist_helper.md)
  : Calculate geographic distance via Haversine formula (helper)
- [`richness()`](https://b-cubed-eu.github.io/dissmapr/reference/richness.md)
  : Calculate Species Richness
- [`turnover()`](https://b-cubed-eu.github.io/dissmapr/reference/turnover.md)
  : Calculate Species Turnover or Beta Diversity
- [`abund()`](https://b-cubed-eu.github.io/dissmapr/reference/abund.md)
  : Calculate Species Abundance
- [`phi_coef()`](https://b-cubed-eu.github.io/dissmapr/reference/phi_coef.md)
  : Calculate Phi Coefficient
- [`cor_spear()`](https://b-cubed-eu.github.io/dissmapr/reference/cor_spear.md)
  : Calculate Spearman's Rank Correlation
- [`cor_pears()`](https://b-cubed-eu.github.io/dissmapr/reference/cor_pears.md)
  : Calculate Pearson's Correlation
- [`diss_bcurt()`](https://b-cubed-eu.github.io/dissmapr/reference/diss_bcurt.md)
  : Calculate Bray-Curtis Dissimilarity
- [`orderwise_diss_gower()`](https://b-cubed-eu.github.io/dissmapr/reference/orderwise_diss_gower.md)
  : Compute Gower dissimilarity between two site vectors
- [`mutual_info()`](https://b-cubed-eu.github.io/dissmapr/reference/mutual_info.md)
  : Calculate Mutual Information
