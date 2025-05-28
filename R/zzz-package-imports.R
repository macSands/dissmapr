# -------------------------------------------------------------------------
#  Internal imports and global-variable declarations
#
#  This file is purely for package hygiene:
#    • collects all @importFrom tags so roxygen can write NAMESPACE cleanly
#    • registers symbols used in data-table / tidy-eval to silence
#      “no visible binding for global variable …” notes.
#
#  It exports **nothing** (hence the terminating NULL).
# -------------------------------------------------------------------------

#' @keywords internal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cor kmeans dist aggregate
#' @importFrom utils combn read.csv unzip
#' @importFrom geosphere distHaversine distm distGeo
#' @importFrom caret findCorrelation
#' @importFrom data.table := rbindlist as.data.table
#' @importFrom dplyr rename mutate select across group_by summarize ungroup
#'             filter relocate bind_cols left_join slice any_of all_of
#'             row_number cur_group_id intersect union
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom pbapply pblapply
#' @importFrom future plan
#' @importFrom cluster daisy pam
#' @importFrom NbClust NbClust
#' @importFrom clValid clValid
#' @importFrom factoextra hcut
#' @importFrom mclust Mclust
#' @importFrom fields Tps
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_fill_gradientn labs
#'             theme_minimal theme element_blank
#' @importFrom sf st_as_sf st_as_sfc st_bbox st_make_grid st_set_geometry
#'             st_centroid st_coordinates st_union st_convex_hull st_buffer
#'             st_join st_within
#' @importFrom terra rast vect ext classify interpolate crds rasterize area
#'             describe intersect union
#' @import     geodata
#' @importFrom geodata worldclim_global
#' @importFrom purrr reduce
#' @importFrom zoo na.approx time<-
#' @importFrom zetadiv Ispline Predict.msgdm Zeta.msgdm
#' @importFrom vegan vegdist
#' @importFrom entropy mi.plugin
#' @importFrom reshape2 melt
NULL

# ---- Silence R CMD check “no visible binding …” --------------------------

utils::globalVariables(c(
  # generic / tidy-eval helpers
  ".", "value", "x", "y", "site_id", "species",

  # data.table assignment symbol
  ":=",

  # temporary or generated column names
  "site_from", "site_to", "pred_zetaExp", "obs_sum",
  "coord_cols", "sp_cols", "orig_x", "orig_y",
  "geometry", "grid_id",

  # mapsheet helpers
  "centroid", "centroid_lat", "centroid_lon",
  "mapsheet", "spp_rich",

  # symbol captured by future::plan inside data.table code
  "plan"
))
