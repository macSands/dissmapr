#’ Internal imports and global variable declarations
#’
#’ This file collects all of our @importFrom directives and silences
#’ “no visible binding” notes via utils::globalVariables().
#’
#’ @keywords internal
#’
#’ @importFrom grDevices    colorRampPalette
#’ @importFrom stats        cor
#’ @importFrom utils        combn read.csv
#’ @importFrom geosphere    distGeo distm
#’ @importFrom caret        findCorrelation
#’ @importFrom terra        centroid area describe intersect union
#’ @importFrom future       plan
#’ @importFrom zetadiv      Ispline Predict.msgdm
#’ @importFrom data.table   `:=` rbindlist as.data.table
#’ @importFrom dplyr        intersect union
#’ @importFrom zoo          time<-
#’
#’ @export
NULL

## avoid “no visible binding” for data.table and other variables:
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "value", "site_from", "site_to", "pred_zetaExp", "obs_sum",
      "coord_cols", "sp_cols", "orig_x", "orig_y", "x", "y",
      ".", "..coord_cols", "..sp_cols", "geometry", "grid_id",
      "map", "centroid", "centroid_lat", "centroid_lon"
    )
  )
}
