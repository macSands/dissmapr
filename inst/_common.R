## vignettes/_common.R   (no YAML, no prose)

# Global chunk defaults -------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.path = "figures/",   # child adds prefix automatically
  message  = FALSE,
  warning  = FALSE,
  echo     = FALSE         # hide code unless a vignette re-enables it
)

# Packages --------------------------------------------------------------
library(dissmapr)   # pkgdown has already installed the package
# library(devtools) # <- *do not* call load_all() inside vignettes
library(httr)
library(geodata)
library(data.table)
library(dplyr)
library(tidyr)
library(zoo)
library(sf)
library(terra)
library(tidyterra)
library(zetadiv)
library(ggplot2)
library(viridis)
library(patchwork)
library(mclust)
library(zen4R)

# Cached objects --------------------------------------------------------
load(
  system.file("extdata", "dissmapr_cache.RData",
              package = "dissmapr"),
  envir = environment()
)
