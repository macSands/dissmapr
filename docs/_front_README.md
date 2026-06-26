# dissmapr

## A Novel Framework for Automated Compositional Dissimilarity and Biodiversity Turnover Analysis

[![repo
status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Release](https://img.shields.io/github/release/b-cubed-eu/dissmapr.svg)](https://github.com/b-cubed-eu/dissmapr/releases)
[![dissmapr status
badge](https://b-cubed-eu.r-universe.dev/dissmapr/badges/version)](https://b-cubed-eu.r-universe.dev/dissmapr)
[![CRAN
status](https://www.r-pkg.org/badges/version/dissmapr)](https://cran.r-project.org//package=dissmapr)
[![R-CMD-check](https://github.com/b-cubed-eu/dissmapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/b-cubed-eu/dissmapr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/b-cubed-eu/dissmapr/graph/badge.svg)](https://app.codecov.io/gh/b-cubed-eu/dissmapr)
[![DOI](https://img.shields.io/badge/DOI-awaiting_upload_to_zenodo-orange)](https://zenodo.org)
[![name status
badge](https://b-cubed-eu.r-universe.dev/badges/:name?color=6CDDB4)](https://b-cubed-eu.r-universe.dev/)
[![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://b-cubed-eu.github.io/dissmapr/LICENSE.md)

------------------------------------------------------------------------

## Introduction

Quantifying **dissimilarity** is central to ecology and biodiversity
science. Whether comparing communities across landscapes, examining
functional distances between species, or assessing environmental
contrasts, dissimilarity metrics underpin analyses of diversity,
turnover, and ecological change. However, these metrics are often
scattered across packages, calculated inconsistently, or tailored to a
single data type (e.g., species composition only).

**`dissmapr`** provides a unified, transparent, and reproducible
framework for calculating and analysing dissimilarities across
**communities**, **traits**, and **environments**. By bringing these
components together, the package helps ecologists explore how
biodiversity varies in space and time, and how functional and
environmental contexts shape these patterns.

The package goes beyond simple pairwise distances, offering tools to:

- Partition dissimilarity into interpretable components such as
  **turnover** and **nestedness**.
- Directly map dissimilarity values to geographic or environmental
  space, enabling spatially explicit analyses.
- Compare metrics across data types, fostering integrated biodiversity
  assessments.

By combining flexible metric choice with consistent data handling and
mapping functions, `dissmapr` is designed to streamline biodiversity
workflows, support cross-study comparability, and provide
decision-relevant outputs for conservation and global change research.

------------------------------------------------------------------------

## Core Concepts

- **Community dissimilarity** - quantifies differences in species
  composition across sites (e.g., Jaccard, Sørensen).
- **Functional dissimilarity** - measures divergence in trait space,
  capturing ecological strategies and redundancy.
- **Environmental dissimilarity** - evaluates contrasts in abiotic
  conditions such as climate, soils, or land cover.
- **Partitioning** - decomposes total dissimilarity into turnover
  (species replacement) vs. nestedness (species loss/gain).
- **Mapping** - links dissimilarity values back to geographic or
  environmental gradients for spatial interpretation.

Together, these concepts provide a flexible toolkit for exploring
biodiversity structure across ecological, functional, and environmental
dimensions.

------------------------------------------------------------------------

## Main Functions

The package is structured around modular steps, so users can build
end-to-end workflows or extract individual components:

### **Data Handling**

- `standardise_traits()` - Clean and harmonise trait datasets for
  consistent use across analyses.
- `compute_distance_matrix()` - Generate dissimilarity matrices from
  traits, environment, or community composition data.

### **Dissimilarity Metrics**

- `beta_diversity()` - Partition dissimilarity into turnover and
  nestedness components, enabling mechanistic interpretation.
- `functional_dissimilarity()` - Calculate distances among species or
  communities in trait space.
- `environmental_dissimilarity()` - Quantify abiotic contrasts between
  sites based on environmental variables.

### **Mapping and Analysis**

- `map_dissimilarity()` - Visualise dissimilarity patterns in geographic
  space, creating reproducible biodiversity maps.
- `summarise_dissimilarity()` - Collapse pairwise dissimilarities into
  site- or region-level summaries.
- `compare_dissimilarity()` - Evaluate agreement or divergence among
  alternative dissimilarity metrics.

------------------------------------------------------------------------

## Applications

`dissmapr` is designed for use across a wide range of biodiversity
research and conservation contexts, including:

- **Community ecology** - comparing species composition across
  landscapes, ecoregions, or time periods.
- **Functional ecology** - linking species trait divergence to
  ecological strategies, redundancy, and resilience.
- **Environmental filtering** - quantifying how abiotic differences
  shape community assembly.
- **Macroecology** - exploring spatial patterns of turnover and
  nestedness at regional to global scales.
- **Conservation planning** - generating reproducible dissimilarity maps
  to identify biodiversity hotspots, refugia, or vulnerable sites.

------------------------------------------------------------------------

## Why Use `dissmapr`?

- **Unified framework** - integrates community, trait, and environmental
  dissimilarities in a consistent manner.
- **Transparent & reproducible** - modular design supports repeatable
  workflows across datasets and studies.
- **Interpretable outputs** - provides partitioned metrics and mapping
  tools for actionable ecological insights.
- **Flexible & extensible** - works with standard R data structures,
  making it easy to integrate into existing pipelines.

------------------------------------------------------------------------

## Installation

The package can be installed from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("b-cubed-eu/dissmapr")
```

------------------------------------------------------------------------

## Citation

If you use **`dissmapr`**, please cite the package and associated
methods. See `citation("dissmapr")` and the repository’s CITATION files.

------------------------------------------------------------------------
