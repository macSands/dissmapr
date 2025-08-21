#' ---
#' title: "Your title"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Your title}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---

#+ setup-0, include=FALSE
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE,
  echo = FALSE, fig.path = "figures/"
)

# Packages --------------------------------------------------------------
library(dissmapr)
library(httr); library(geodata); library(data.table); library(dplyr)
library(tidyr); library(zoo); library(sf); library(terra); library(tidyterra)
library(zetadiv); library(ggplot2); library(viridis); library(patchwork)
# guard optional deps
has_zen4r <- requireNamespace("zen4R", quietly = TRUE)
if (has_zen4r) library(zen4R)
library(purrr); library(RColorBrewer)

#+ cache, include=FALSE
load(system.file("extdata", "dissmapr_cache.RData", package = "dissmapr"),
     envir = knitr::knit_global())
